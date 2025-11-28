import sys; sys.path.append('../3rdparty/ElasticRods/python')
import elastic_rods, elastic_knots
import numpy as np, matplotlib.pyplot as plt, time, io, os

from helpers import *
from parametric_curves import *
import py_newton_optimizer

from linkage_vis import LinkageViewer as Viewer, CenterlineViewer
from tri_mesh_viewer import PointCloudViewer, PointCloudMesh

import parallelism
parallelism.set_max_num_tbb_threads(1)

class MEP:
    def __init__(
            self,
            start_file,
            goal_file,
            size_reduction = 4,
            rod_radius = 0.2,
            material_stiffness = 2000,
            poisson_ratio = 0.3,
            contact_stiffness=10000,
            hasCollisions = True
            
    ):
        """
        Initiate Minimum Energy Path class with all necessary parameters.
        Args:
            start_file (string): Path to the start file.
            goal_file (string): Path to the goal file.
            size_reduction (float): Size reduction factor for the knot. Defaults to 4.
            rod_radius (float): Radius of the rod. Defaults to 0.2.
            material_stiffness (float): Stiffness of the material. Defaults to 2000.
            poisson_ratio (float): Poisson ratio. Defaults to 0.3.
            contact_stiffness (float): Stiffness of the contact. Defaults to 10000.
        Notes:
              - the default parameters represent a steel wire with 0.2mm radius
        """
        self.start_file = start_file
        self.goal_file = goal_file

        self.rod_radius = rod_radius
        material = elastic_rods.RodMaterial('ellipse', material_stiffness, poisson_ratio, [rod_radius, rod_radius])
        self.start = read_nodes_from_file(start_file)[::size_reduction]
        self.goal = read_nodes_from_file(goal_file)[::size_reduction]
        pr = define_periodic_rod(self.start, material)
        self.rod_list = elastic_knots.PeriodicRodList([pr])
        self.nPoints = len(self.start)

        self.optimizerOptions = py_newton_optimizer.NewtonOptimizerOptions()
        self.optimizerOptions.niter = 10000
        self.optimizerOptions.gradTol = 1e-8
        self.hessianShift = 1e-4 * compute_min_eigenval_straight_rod(pr)

        self.problemOptions = elastic_knots.ContactProblemOptions()
        self.problemOptions.hasCollisions = hasCollisions
        self.problemOptions.contactStiffness = contact_stiffness
        self.problemOptions.dHat = 2 * rod_radius

        self.contactProblem = elastic_knots.ContactProblem(self.rod_list, self.problemOptions)

        self.path = []

    def loadPath(self,file):
        loaded = np.load(file)
        path = [loaded[key] for key in loaded.files]
        return path

    def reducePath(self,path,factor):
        return path[::factor] + [path[-1]]

    def savePath(self,path,file):
        np.savez(file, *path)

    def getPath(self):
        return self.path

    def setPath(self, path):
        self.path = path

    def showPath(self,path,wait=0.05):
        for DoFs in path:
            self.rod_list.setDoFs(DoFs)
            time.sleep(wait)
            self.view.update()

    def getPathEnergy(self,path):
        path_energy = []
        for DoFs in path:
            self.contactProblem.setDoFs(DoFs)
            path_energy.append(self.contactProblem.energy())
        return path_energy

    def getContactProblem(self):
        return self.contactProblem

    def getRodList(self):
        return self.rod_list

    def setViewer(self, viewer):
        """
        set a Viewer. The ElasticKno viewer only works in jupyter Notebooks.
        """
        self.view = viewer

    def unknottingByRandomizing(
            self,
            iterations = 50000,
            it_rand = 1500,
            it_target = 2500,
            holding_dist = 4.0,
            solution_dist = 100,
            move_size = 0.1,
            step_size = 0.01
    ):
        """
        Attempts to reach the goal configuration by randomizing intermediate targets.

        Args:
            iterations (int, optional): Total number of iterations to run. Defaults to 50,000.
            it_rand (int, optional): Iteration at which to switch to a random target. Defaults to 1,500.
            it_target (int, optional): Iteration at which to return to the original goal. Defaults to 2,500.
            holding_dist (float, optional): Distance threshold to "hold" a point during randomization. Defaults to 4.0.
            solution_dist (float, optional): Distance threshold to consider the goal reached. Defaults to 100.0.
            move_size (float, optional): Scaling factor for the movement force toward the target. Defaults to 0.1.
            step_size (float, optional): Step size for updating positions. Defaults to 0.01.

        Returns:
            self.path
        Notes:
            - Assumes self.goal, self.nPoints, self.contactProblem, and self.view are defined.
            - Stores intermediate configurations in `self.path`.
        """

        def moveTo(p, target):
            F = 2 * (target - p)
            return np.array(F)

        def getSpringForce(R, R_next, k=2.0):
            d = R_next - R
            F = k * d
            return F

        #not used, but cool concept
        # when there is contact slide the edges along each other
        def sliding(a, b, target):
            o = np.cross(a, b)
            if np.dot(o, target) < 0:
                o = -o
            return o / np.linalg.norm(o)

        #try to find a path by randomizing inbetween
        j = 0
        self.path = []
        start_DoFs = self.contactProblem.getDoFs().copy()
        start_DoFs[:self.nPoints *3]=self.start.flatten()
        self.contactProblem.setDoFs(start_DoFs)

        for it in range(iterations):
            DoFs = self.contactProblem.getDoFs().copy()
            points = DoFs[:self.nPoints * 3].reshape((-1, 3))
            grad = self.contactProblem.contactForces()
            contactForce = grad[:self.nPoints * 3].reshape((-1, 3))

            if j == 0:
                target = self.goal

            elif j == it_rand:
                minVal = np.min(self.goal)
                maxVal = np.max(self.goal)
                randVal = np.random.rand(*self.goal.shape) * (maxVal - minVal) + minVal
                for i in range(self.nPoints):
                    if np.linalg.norm(points[i] - target[i]) < holding_dist:
                        randVal[i] = target[i]
                        print(f"holding {i}")

                target = randVal
            elif j == it_target:
                j = -1
            j += 1

            for i in range(0, self.nPoints):

                if i == self.nPoints - 1:
                    springForce = getSpringForce(points[i], points[0]) + getSpringForce(points[i], points[i - 1])
                else:
                    springForce = getSpringForce(points[i], points[i + 1]) + getSpringForce(points[i], points[i - 1])

                F = move_size * moveTo(points[i], target[i]) + contactForce[i] + springForce

                F = F / np.linalg.norm(F)
                points[i] = (points[i] + step_size * F)

            DoFs[:self.nPoints * 3] = points.flatten()
            self.contactProblem.setDoFs(DoFs.flatten())

            if it % 100 == 0:
                print(f"{it} Dist to Target: {np.linalg.norm(points - target)}")

            if it % 25 == 0:
                self.path.append(self.contactProblem.getDoFs())

            if (np.linalg.norm(points - self.goal) < solution_dist):
                print("Found a solution")
                break

            self.view.update()
        return self.path

    def relaxPathEnd(self, path):
        """
        relaxes the end of a given path to reach the goal configuration.
        Args:
            path (list, optional): List of paths to relax. Defaults to self.path.
        Returns:
            totalPath (list): The given Path + the relaxation
        Notes:
            - if you like the result, use setPath(totalPath)
        """
        self.contactProblem.setDoFs(path[-1])
        relaxPath = []
        self.view.update()
        def callback(problem, iteration):
            if iteration % 20 == 0:
                relaxPath.append(self.rod_list.getDoFs())
                self.view.update()

        report = elastic_knots.compute_equilibrium(
            self.rod_list, self.problemOptions, self.optimizerOptions,
            externalForces=np.zeros(self.rod_list.numDoF()),
            softConstraints=[],
            callback=callback,
            hessianShift=self.hessianShift
        )
        view.update()

        totalPath = path + relaxPath

        return totalPath

    def findSaddlePoints(
            self,
            path,
            i_start = 0,
            i_goal = -1,
            wait = 1
    ):
        """
        This funktion helps you to identify saddle points by relaxing every knot in the path.
        Identifying is done by hand.
        Args:
            i_start (int, optional): Starting index. Defaults to 0.
            i_goal (int, optional): Goal index. Defaults to -1.
            path (list, optional): The path to identify saddle points within. Defaults to self.path.
        Returns:
            shortend path.
        Notes:
            - it helps to shorten the path
            - use i_start and i_goal to narrow the search down
            - in the end i_start should be the last knot which relaxes to the start state
            - in the end i_goal should be the first knot which relaxes to the goal state

        """
        shortPath = path[i_start:i_goal]
        for i, Dof in zip(range(len(shortPath)), shortPath):
            relaxPath =[]
            self.contactProblem.setDoFs(Dof)

            def callback(problem, iteration):
                if iteration % 10 == 0:
                    relaxPath.append(self.rod_list.getDoFs())
                    self.view.update()

            report = elastic_knots.compute_equilibrium(
                self.rod_list, self.problemOptions, self.optimizerOptions,
                externalForces=np.zeros(self.rod_list.numDoF()),
                softConstraints=[],
                callback=callback,
                hessianShift=self.hessianShift
            )
            if i == 0:
                reversePath = relaxPath[-1::-1]
            self.view.update()
            time.sleep(wait)
        return reversePath[::4] + shortPath + relaxPath[::4] + [relaxPath[-1]]

    def optimizePath(
            self,
            R,
            step_size=0.00002,
            min_step = 1e-8,
            max_step = 1e-3,
            iterations=10000,
            gradTol = 3e-3
    ):
        """
        Optimize the given path using the NEB force given by TODO link
        Args:
            R (list, optional): List of paths to optimize. Defaults to self.path.
            rod (list, optional): the contactProblem. Defaults to self.contactProblem.
            step_size (float, optional): Step size for updating positions. Defaults to 0.00002.
            min_step (float, optional): Minimum step size for updating positions. Defaults to 1e-5.
            max_step (float, optional): Maximum step size for updating positions. Defaults to 1e-2.
            iterations (int, optional): Number of iterations. Defaults to 10000.
            gradTol (float, optional): Tolerance for gradient norms. Defaults to 3e-3.
        Returns:
            R (list): optimized path
        """

        def safe_normalize(v, eps=1e-12):
            norm = np.linalg.norm(v)
            if norm < eps or np.isnan(norm):
                return np.zeros_like(v)
            return v / norm

        def getTangent(R_pre, R_next):
            t = R_next - R_pre
            t = safe_normalize(t)
            return np.array(t)

        # F_spring = k * (||R_{i+1} - R_i|| - ||R_i - R_{i-1}||) * tangent_i
        def getSpringForce(R_pre, R, R_next, k = 2):
            tangent = getTangent(R_pre, R_next)
            F = k * (np.linalg.norm(R_next - R) - np.linalg.norm(R - R_pre)) * tangent
            return np.array(F)

        # F_perp =-d_R_{i} + d_R_{i}°tangent_i*tanget_i
        def getPerpForce(d_R, R_pre, R_next):
            tangent = getTangent(R_pre, R_next)
            F = (-d_R + (np.dot(d_R, tangent) * tangent))
            return np.array(F)

        def getNEBForce(d_R, R_pre, R, R_next):
            return getPerpForce(d_R, R_pre, R_next) + getSpringForce(R_pre, R, R_next)

        # Initialize per-knot step sizes
        step_size = np.full(len(R), step_size)
        rod = self.contactProblem
        prevE = np.full(len(R), np.inf)
        start_time = time.time()
        for it in range(iterations + 1):
            totalE = 0.0
            totalGradNorm = 0
            for i in range(1, len(R) - 1):  # skip fixed endpoints
                rod.setDoFs(R[i])
                d_R = rod.gradient()
                F = getNEBForce(d_R, R[i - 1], R[i], R[i + 1])
                E = rod.energy()
                totalE += E
                totalGradNorm += np.linalg.norm(d_R)
                # Compare force magnitudes to adapt step size per knot
                if E < prevE[i]:
                    step_size[i] *= 1.05  # Gradually increase
                else:
                    step_size[i] *= 0.5  # Back off if divergence

                step_size[i] = np.clip(step_size[i], min_step, max_step)

                prevE[i] = E

                R[i] += step_size[i] * F  # Gradient ascent — be sure F is negative grad if you want descent

            if (totalGradNorm / (len(R) - 2)) < gradTol:
                print(f"Converged at iteration {it}")
                break

            if it % 100 == 0:
                print(f"{it} totalEnergy = {totalE:.6f}, after {(time.time() - start_time):.2f}s")
                print(
                    f"  step_size (min/max): {step_size.min():.2e} / {step_size.max():.2e} mean: {np.mean(step_size)}")
                print(f"  total gradient normalised: {totalGradNorm}")

        return R





