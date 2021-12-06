import math

from django.contrib.gis.geos import Point
import numpy as np


class RegionOfSupport:
    """
    Class based on the implementation of RegionOfSupport.java by Luı́s Carlos in his thesis "Morphing techniques in
    spatio-temporal database" ("Aplicação de técnicas de morphing em bases de dados espacio-temporais")
    """

    def __init__(self, points, fp):
        self.center = Point(0, 0)
        self.points = points                # polygon points
        self.covariance = np.zeros((2, 2))
        self.e_N, self.ev_N = None, 0       # Normal Vector and Normal Value
        self.e_T, self.ev_T = None, 0       # Tangent Vector and Tangent Value

        if len(points) > 0:
            self.h = len(points) / 2
            self.fp = fp                    # feature point
            self.calculate_center()
            self.calculate_covariance()
            e = self.calculate_eigen_values_and_eigen_vectors()
            self.calculate_normal_and_tangent_vectors(
                [e[1][0], e[2][0]], e[0][0],
                [e[1][1], e[2][1]], e[0][1]
            )
        else:
            print("Error in RoS: there are no points")

    def calculate_center(self):
        x, y = 0, 0
        for point in self.points:
            if isinstance(point, tuple):
                x += point[0]
                y += point[1]
            else:
                x += point.point[0]
                y += point.point[1]
        x = x / (2 * self.h + 1)
        y = y / (2 * self.h + 1)
        self.center = Point(x, y)

    def calculate_covariance(self):
        for point in self.points:
            if isinstance(point, tuple):
                x = point[0]
                y = point[1]
            else:       # in the case of being a feature point
                x = point.point[0]
                y = point.point[1]
            x -= self.center[0]
            y -= self.center[1]

            self.covariance[0][0] += x * x
            self.covariance[1][0] = x * y
            self.covariance[0][1] += x * y
            self.covariance[1][1] += y * y

        self.covariance[0][0] /= 2 * self.h + 1
        self.covariance[0][1] /= 2 * self.h + 1
        self.covariance[1][0] /= 2 * self.h + 1
        self.covariance[1][1] /= 2 * self.h + 1

    def calculate_eigen_values_and_eigen_vectors(self):
        ev1, ev2, res = [], [], np.zeros((3, 2))
        a, b, c, d = self.covariance[0][0], self.covariance[0][1], self.covariance[1][0], self.covariance[1][1]
        T = a + d
        D = a * d - b * c
        # TODO: check this (abs not in the original code)
        L1 = T / 2 + math.sqrt(abs((T * T) / 4 - D))
        L2 = T / 2 - math.sqrt(abs((T * T) / 4 - D))

        if c != 0:
            ev1.append(L1 - d)
            ev1.append(c)

            ev2.append(L2 - d)
            ev2.append(c)
        elif b != 0:
            ev1.append(b)
            ev1.append(L1 - a)

            ev2.append(b)
            ev2.append(L2 - a)
        else:
            ev1.append(1)
            ev1.append(0)

            ev2.append(0)
            ev2.append(1)

        # return the results
        np.zeros((3, 2))
        res[0][0], res[1][0], res[2][0] = L1, ev1[0], ev1[1]
        res[0][1], res[1][1], res[2][1] = L2, ev2[0], ev2[1]
        return res

    def calculate_normal_and_tangent_vectors(self, eigen_vector1, eigen_value1, eigen_vector2, eigen_value2):
        bisector = self.calculate_bisector()

        # Calculate the dot product for each EigenVector
        # dot product <v1, v2> = v1,1*v2,1 + v1,2*v2,2 + ... + v1,n*v2,n
        # Center the vectors in FP

        nX = eigen_vector1[0] + (self.fp.point[0] - self.center[0])
        nY = eigen_vector1[1] + (self.fp.point[1] - self.center[1])
        dot_ev1 = nX * bisector[0] + nY * bisector[1]

        nX = eigen_vector2[0] + (self.fp.point[0] - self.center[0])
        nY = eigen_vector2[1] + (self.fp.point[1] - self.center[1])
        dot_ev2 = nX * bisector[0] + nY * bisector[1]

        # choose e_T and e_N
        e_N, e_T = np.zeros((2, 2)),  np.zeros((2, 2))
        e_N[0][0], e_N[0][1] = self.center[0], self.center[1]
        e_T[0][0], e_N[0][1] = self.center[0], self.center[1]

        # the eigenvector with the smallest dot product is considered the tangent vector
        if dot_ev1 > dot_ev2:
            # eigen_vector1 is the normal vector
            e_N[1][0], e_N[1][1] = self.center[0] + eigen_vector1[0], self.center[1] + eigen_vector1[1]
            ev_N = eigen_value1

            e_T[1][0], e_T[1][1] = self.center[0] + eigen_vector2[0], self.center[1] + eigen_vector2[1]
            ev_T = eigen_value2
        else:
            # eigen_vector2 is the normal vector
            e_N[1][0], e_N[1][1] = self.center[0] + eigen_vector2[0], self.center[1] + eigen_vector2[1]
            ev_T = eigen_value2

            e_T[1][0], e_T[1][1] = self.center[0] + eigen_vector1[0], self.center[1] + eigen_vector1[1]
            ev_N = eigen_value1

        self.e_N, self.ev_N = e_N, ev_N
        self.e_T, self.ev_T = e_T, ev_T

    def calculate_bisector(self):
        _prev, _, _next = self.fp.get_least_opening_triangle_points()
        bisector = [_prev[0] + abs(_prev[0] - _next[0]) / 2, _prev[1] + abs(_prev[1] - _next[1]) / 2]
        bisector[0] -= self.fp.point[0]
        bisector[1] -= self.fp.point[1]
        return bisector




