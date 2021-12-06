import math
import numpy as np
from correspondence import Correspondence


class PolygonCorrespondence:
    """
    Class based on the implementation of PolygonCorrespondences.java by Luı́s Carlos in his thesis "Morphing techniques
    in spatio-temporal database" ("Aplicação de técnicas de morphing em bases de dados espacio-temporais")
    """

    def __init__(self, s_fps, t_fps, ft_variation, ft_side, ft_size, weight_sim):
        self.s_fps = s_fps  # source polygon feature points
        self.t_fps = t_fps  # target polygon feature points
        self.weight_sim = weight_sim  # weight of the similarity value of two feature points

        # calculate the similarity costs
        self.sim_costs = np.zeros((len(s_fps), len(t_fps)))
        for s in range(len(s_fps)):
            for t in range(len(t_fps)):
                self.sim_costs[s][t] = s_fps[s].calculate_similarity(t_fps[t], ft_variation, ft_side, ft_size)

        # calculate the discard costs
        self.discard_cost_s, self.discard_cost_t = [], []
        for i in range(len(s_fps)):
            self.discard_cost_s.append(s_fps[i].calculate_discard_cost(ft_variation, ft_side, ft_size))
        for j in range(len(t_fps)):
            self.discard_cost_t.append(t_fps[j].calculate_discard_cost(ft_variation, ft_side, ft_size))

    def get_feature_point_correspondences(self, skips):
        _min = math.inf
        correspondences = []

        for s_init in range(len(self.s_fps)):
            for t_init in range(len(self.t_fps)):
                tmp = self.get_feature_point_correspondences_by_path(s_init, t_init, skips)
                total = 0
                for c in tmp:
                    total += c.cost

                if _min > total:
                    _min = total
                    correspondences = tmp
        return correspondences

    def get_path(self, s_init, t_init, skips):
        correspondences = []
        s_i, t_i = s_init, t_init
        delta_s, delta_t = 0, 0

        while True:
            # direct correspondence
            s_next, t_next = s_i, t_i
            _min = self.delta(s_i, s_next, t_i, t_next)

            # skips
            nr_s_skips, nr_t_skips = 0, 0  # number os skips
            for s_skips in range(skips + 1):
                for t_skips in range(skips + 1):
                    tmp = self.delta(s_i, (s_i + s_skips) % len(self.s_fps), t_i, (t_i + t_skips) % len(self.t_fps))
                    if tmp < _min:
                        _min = tmp
                        nr_s_skips = s_skips
                        nr_t_skips = t_skips
                        s_next, t_next = (s_i + s_skips) % len(self.s_fps), (t_i + t_skips) % len(self.t_fps)

            # increment already processed vertexes
            delta_s += nr_s_skips
            delta_t += nr_t_skips

            if delta_s > len(self.s_fps) or delta_t > len(self.t_fps):
                break

            # Add correspondence
            c = Correspondence()
            c.point_s = self.s_fps[s_next]
            c.s_i = s_next
            c.point_t = self.t_fps[t_next]
            c.t_i = t_next
            c.cost = _min
            correspondences.append(c)

            # move to next values
            s_i, t_i = (s_next + 1) % len(self.s_fps), (t_next + 1) % len(self.t_fps)
            delta_s += 1
            delta_t += 1

            if delta_s >= len(self.s_fps) or delta_t >= len(self.t_fps):
                break

        return correspondences

    def delta(self, s_begin, s_end, t_begin, t_end):
        s_begin, s_end = s_begin % len(self.s_fps), s_end % len(self.s_fps)
        t_begin, t_end = t_begin % len(self.t_fps), t_end % len(self.t_fps)
        s_dist, t_dist = abs(s_end - s_begin), abs(t_end - t_begin)

        costs = 0
        for i in range(s_dist):
            costs += self.discard_cost_s[(s_begin + i) % len(self.s_fps)]
        for j in range(t_dist):
            costs += self.discard_cost_t[(t_begin + j) % len(self.t_fps)]
        costs += self.weight_sim * self.sim_costs[s_end][t_end]
        return costs

    def get_feature_point_correspondences_by_path(self, s_init, t_init, skips):
        return self.get_path(s_init, t_init, skips)
