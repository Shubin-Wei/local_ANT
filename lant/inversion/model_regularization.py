"""lant's regularization of the 3-D inversion (inversion.regularization: lant), in place of DSurfTomo's.

Each iteration still solves for the update dv of the current model m (rays and kernels in m, as DSurfTomo), but the
penalties act on the model itself, m + dv - m0 (m0 the starting model), not on the update:

    min  (1/N) ||W (G dv - r)||^2  +  (w^2/n) ||L (m + dv - m0)||^2

- "jumping" instead of "creeping" (Shaw & Orcutt 1985, Geophys. J. R. astr. Soc. 82, 375-414): the smoothing holds
  the final model, so the weight sets how rough the result may be, whatever the number of iterations; DSurfTomo
  damps and smooths each update, whose sum is not controlled.  The same choice as Occam's
  inversion (Constable, Parker & Constable 1987, Geophysics 52, 289-300) and as lant's 1-D inversion.
- L: the 7-point Laplacian of DSurfTomo's interior rows (6 at the node, -1 at each neighbour, grid-index units: the
  lateral spacing is equal in km in both directions, the depth nodes follow the resolution), extended to the nodes on
  the faces of the parameter box with the neighbours they have (the graph Laplacian; a one-sided difference at a
  face; Aster, Borchers & Thurber 2013, Parameter Estimation and Inverse Problems, 2nd ed., ch. 4 roughening
  matrices).  DSurfTomo gives the face nodes, the top layer included, only 2w on the diagonal and no smoothing.
- 1/N: the data term is the mean squared residual of the N rows kept by the rejection rule (a misfit per datum,
  Aster et al. 2013 ch. 2; Constable et al. 1987), so the weight does not change with the number of travel times.
  With DSurfTomo's sum over all rows the same weight is weaker the more travel times there are.  For the same
  reason the model terms are means over the n parameters: w is then the ratio of the rms travel-time misfit (s) to
  the rms roughness (km/s) it trades against, the same number for any grid and data set.

- damping (the damp of DSurfTomo.in): rows d (m + dv - m0) / sqrt(n), towards the starting model, also as a mean
  over the parameters.  Nodes no ray reaches are held by it near the starting model instead of following the
  smooth continuation of their neighbours without any data (zeroth-order Tikhonov term; Aster et al. 2013 ch. 4).

    min  (1/N) ||W (G dv - r)||^2  +  (w^2/n) ||L (m + dv - m0)||^2  +  (d^2/n) ||m + dv - m0||^2

The data rows (G, rejection) and the solver (LSMR, without its own damping) are DSurfTomo's.
"""

import numpy as np

F = np.float32


def laplacian(nvx, nvz, nzu):
    """Graph Laplacian of the parameter box (nzu depths x nvz longitudes x nvx latitudes, node index
    k nvz nvx + j nvx + i as DSurfTomo's columns) -> (rows, cols, values) with unit edge weights."""
    k, j, i = np.meshgrid(np.arange(nzu), np.arange(nvz), np.arange(nvx), indexing="ij")
    k, j, i = k.ravel(), j.ravel(), i.ravel()
    node = k * nvz * nvx + j * nvx + i
    rows, cols, vals = [node], [node], [np.zeros(node.size, F)]
    for has, step in ((i > 0, -1), (i < nvx - 1, 1), (j > 0, -nvx), (j < nvz - 1, nvx),
                      (k > 0, -nvz * nvx), (k < nzu - 1, nvz * nvx)):  # fmt: skip
        rows.append(node[has])
        cols.append(node[has] + step)
        vals.append(np.full(int(has.sum()), -1.0, F))
        vals[0][has] += F(1.0)
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    order = np.lexsort((cols, rows))  # grouped by row, as the LSMR matrix stores them
    return rows[order], cols[order], vals[order]


class Regularization:
    """The fixed parts of lant's regularization of one project: L, the starting model and the weight."""

    def __init__(self, nvx, nvz, nzu, m0, weight, damping=0.0):
        self.n = nvx * nvz * nzu
        self.l_rows, self.l_cols, self.l_vals = laplacian(nvx, nvz, nzu)
        self.m0 = np.asarray(m0, np.float64).ravel()
        self.weight = float(weight)
        self.damping = float(damping)
        self.per_node = 1.0 / np.sqrt(self.n)  # the model terms as means over the n parameters

    def rows(self, m, first_row):
        """Smoothing (and damping) rows below the data rows: (rows, cols, values, rhs) for the update of model m."""
        x = np.asarray(m, np.float64).ravel() - self.m0
        lx = np.zeros(self.n)
        np.add.at(lx, self.l_rows, self.l_vals.astype(np.float64) * x[self.l_cols])
        w = self.weight * self.per_node
        rows, cols = first_row + self.l_rows, self.l_cols
        vals, rhs = (self.l_vals * F(w)).astype(F), (-w * lx).astype(F)
        if self.damping > 0:
            d = self.damping * self.per_node
            node = np.arange(self.n)
            rows = np.r_[rows, first_row + self.n + node]
            cols = np.r_[cols, node]
            vals = np.r_[vals, np.full(self.n, d, F)]
            rhs = np.r_[rhs, (-d * x).astype(F)]
        return rows, cols, vals, rhs

    def row_count(self):
        return 2 * self.n if self.damping > 0 else self.n
