import collinearw
from collinearw import Histogram, Histogram2D, HistMaker
from collinearw.core import HistogramBase
from tqdm import tqdm
import awkward as ak
import numpy as np
import logging
import copy
import numbers
from awkward._connect import numexpr


log = logging.getLogger(__name__)

ne_evaluate = numexpr.evaluate


class Graph(HistogramBase):
    def __init__(self, name, x, y, xtitle, ytitle, type, filter_type=None):
        super().__init__(name)
        self.xvar = x
        self.yvar = y
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.xdata = []
        self.ydata = []
        self.limit = 20
        self.counter = 0
        self.type = type
        self.reach_limit = False
        self.filter_type = filter_type
        self.action_before_fill = None

    def __copy__(self):
        c_self = super().__copy__()
        c_self.xvar = self.xvar
        c_self.yvar = self.yvar
        c_self.xtitle = self.xtitle
        c_self.ytitle = self.ytitle
        c_self.limit = self.limit
        c_self.counter = self.counter
        c_self.type = self.type
        c_self.reach_limit = self.reach_limit
        c_self.filter_type = self.filter_type
        c_self.action_before_fill = self.action_before_fill
        c_self.xdata = copy.copy(self.xdata)
        c_self.ydata = copy.copy(self.ydata)
        return c_self

    def __deepcopy__(self, memo):
        c_self = super().__deepcopy__(memo)
        c_self.xvar = self.xvar
        c_self.yvar = self.yvar
        c_self.xtitle = self.xtitle
        c_self.ytitle = self.ytitle
        c_self.limit = self.limit
        c_self.counter = self.counter
        c_self.type = self.type
        c_self.reach_limit = self.reach_limit
        c_self.filter_type = self.filter_type
        c_self.action_before_fill = self.action_before_fill
        c_self.xdata = copy.deepcopy(self.xdata, memo)
        c_self.ydata = copy.deepcopy(self.ydata, memo)
        return c_self

    def copy(self, *args, **kwargs):
        return copy.deepcopy(self)

    def add(self, rhs):
        self.xdata = self.xdata + rhs.xdata
        self.ydata = self.ydata + rhs.ydata

    def __add__(self, rhs):
        c_self = copy.deepcopy(self)
        c_self.add(rhs)
        return c_self

    @property
    def hist_type(self):
        return "graph"

    @property
    def observable(self):
        return (self.xvar, self.yvar)

    @property
    def ndata(self):
        return len(self.xdata)

    def from_array(self, xdata, ydata):
        if self.counter == self.limit:
            self.reach_limit = True
            return
        if self.action_before_fill == "flatten":
            self.xdata.append(ak.flatten(xdata).to_numpy())
            self.ydata.append(ak.flatten(ydata).to_numpy())
            self.counter += 1
            return
        for x, y in zip(xdata, ydata):
            self.xdata.append(x.to_numpy())
            self.ydata.append(y.to_numpy())
            self.counter += 1


class SSRLHisto1D(Histogram):
    def __init__(self, *args, selection=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.selection = collinearw.utils.to_numexpr(selection)

    def __copy__(self):
        c_self = super().__copy__()
        c_self.selection = self.selection
        return c_self

    def __deepcopy__(self, memo):
        c_self = super().__deepcopy__(memo)
        c_self.selection = copy.deepcopy(self.selection, memo)
        return c_self

    def from_array(self, event, mask=None, func=None):
        if func is None:
            if self.selection is not None:
                hist_mask = ne_evaluate(self.selection, event)
                breakpoint()
        else:
            raise NotImplementedError
            data = func(event, self.selection)
        # super().from_array(data)


class AvgGraph(Graph):

    @property
    def hist_type(self):
        return "avg-graph"

    def from_array(self, xdata, ydata):
        if self.reach_limit:
            return
        if self.counter == self.limit:
            self.reach_limit = True
            return
        if self.action_before_fill == "flatten":
            self.xdata.append(ak.flatten(xdata).to_numpy())
            self.ydata.append(ak.flatten(ydata).to_numpy())
            self.counter += 1
            return
        if self.xdata == []:
            # self.xdata = np.average(xdata.to_numpy(), axis=0)
            self.xdata = np.sum(xdata.to_numpy(), axis=0)
        else:
            self.xdata += np.sum(xdata.to_numpy(), axis=0)
        if self.ydata == []:
            self.ydata = np.sum(ydata.to_numpy(), axis=0)
        else:
            self.ydata += np.sum(ydata.to_numpy(), axis=0)
        self.counter += len(xdata)


class SSRLHistMaker(HistMaker):

    early_termination_counter = 20

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.histogram_1d_alias = []

    def plevel_process(self, p, file_name, *, branch_list=None):
        with self.open_file(file_name) as tfile:
            ttree = tfile[p.treename]
            # check ttree and number of entry. return early if it's zero
            if ttree is None or ttree.num_entries == 0:
                return p

            if p.weights:
                if isinstance(p.weights, list):
                    process_weights = "*".join(p.weights)
                else:
                    process_weights = p.weights
            else:
                process_weights = None
            log.debug(f"Process level weights: {process_weights}")

            # try to get branches from the Process instance, and the branches of
            # regions within it.
            p_branch = p.ntuple_branches.copy()
            for r in p.regions:
                p_branch |= r.ntuple_branches
            if self.branch_list is not None:
                p_branch |= set(self.branch_list)

            branch_filter = branch_list or p_branch or self.branch_list
            # if branch renaming is requested, we need to make sure the original
            # names are used for branch filtering.
            if self.branch_rename:
                for old_bname, new_bname in self.branch_rename.items():
                    if new_bname not in branch_filter:
                        continue
                    branch_filter.discard(new_bname)
                    branch_filter.add(old_bname)

            with tqdm(
                desc=f"Processing {p.name}",
                total=ttree.num_entries,
                leave=False,
                unit="events",
                dynamic_ncols=True,
                disable=self.disable_pbar,
            ) as pbar_events:
                for event, report in ttree.iterate(
                    step_size=self.step_size,
                    filter_name=branch_filter,
                    report=True,
                    # library="np",
                ):
                    nevent = report.tree_entry_stop - report.tree_entry_start
                    pbar_events.set_description(f"Processing {nevent} events")

                    # all_mask is a mask with only process level selection
                    # if no process level seletion, accept all events.
                    # the mask is array of True/False.
                    if p.selection_numexpr:
                        all_mask = ne_evaluate(p.selection_numexpr, event)
                        if not ak.count_nonzero(all_mask):
                            log.debug("No event after process selection")
                            pbar_events.update(nevent)
                            continue
                    else:
                        all_mask = None  # ak.ones_like(event)

                    pbar_regions = tqdm(
                        p.regions,
                        leave=False,
                        unit="regions",
                        disable=self.disable_pbar,
                    )
                    for r in pbar_regions:
                        pbar_regions.set_description(
                            f"{p.name}, Region: {r.name}({len(r.histograms)})"
                        )

                        # setting process level and region level selection
                        # this is basically cuts in ROOT TTree but in numpy format
                        p_r_selection = []
                        if p.selection_numexpr:
                            p_r_selection.append(f"({p.selection_numexpr})")
                        if r.selection_numexpr:
                            p_r_selection.append(f"({r.selection_numexpr})")

                        # combinding process and region level selections
                        if p_r_selection:
                            selection_str = "&".join(p_r_selection)
                            selection_str = selection_str.replace("()", "")
                            selection_str = selection_str.strip().strip("&")
                        else:
                            selection_str = ""

                        # if no selection string is found, assume accepting all values
                        if selection_str:
                            mask = ne_evaluate(selection_str, event)
                        else:
                            log.debug(
                                f"empty seleciton on region {r.name}. Assume no selection."
                            )
                            mask = all_mask

                        # obtaining process and region level weights
                        # if none of them were found, try to use the
                        # histmaker default weight. i.e self.default_weight
                        weights = None
                        if process_weights is None and r.weights is None:
                            if self.default_weight:
                                weights = ne_evaluate(self.default_weight, event)
                        else:
                            if process_weights:
                                weights = ne_evaluate(process_weights, event)
                            if r.weights:
                                if isinstance(r.weights, list):
                                    for w in r.weights:
                                        if weights is None:
                                            weights = event[w]
                                        else:
                                            weights *= event[w]
                                elif str.isnumeric(r.weights):
                                    if weights is None:
                                        weights = float(r.weights)
                                    else:
                                        weights *= float(r.weights)
                                else:
                                    mul_w = ne_evaluate(r.weights, event)
                                    if weights is None:
                                        weights = mul_w
                                    else:
                                        weights *= mul_w
                            # check if user enforce to use default weight
                            if self.enforce_default_weight and self.default_weight:
                                weights *= ne_evaluate(self.default_weight, event)

                        if mask is None:
                            m_mask = None
                        elif mask.ndim == 1:
                            m_mask = mask
                        else:
                            m_mask = ak.flatten(mask)

                        if weights is None:
                            w = None
                        elif isinstance(weights, numbers.Number):
                            if m_mask is None:
                                w = weights * mask
                            else:
                                w = weights * m_mask
                        else:
                            if m_mask is None:
                                w = ak.flatten(weights)[m_mask]
                            else:
                                w = ak.flatten(weights)

                        for hist in r.histograms:
                            if hist.hist_type == '2d':
                                xobs, yobs = hist.observable
                                xdata = ne_evaluate(xobs, event)
                                ydata = ne_evaluate(yobs, event)
                                histw = None
                                if mask is not None:
                                    xdata = xdata[mask]
                                    ydata = ydata[mask]
                                    # if mask is not None:
                                    #     xdata = xdata[ak.any(mask,axis=1)]
                                    #     ydata = ydata[ak.any(mask,axis=1)]
                                    if w is not None:
                                        histw = w[m_mask]
                                if xdata.ndim != 1:
                                    xdata = ak.flatten(xdata)
                                if ydata.ndim != 1:
                                    ydata = ak.flatten(ydata)
                                if ak.any(xdata) and ak.any(ydata):
                                    hist.from_array(xdata, ydata, histw)
                            elif hist.hist_type == "graph":
                                xobs, yobs = hist.observable
                                xdata = ne_evaluate(xobs, event)
                                ydata = ne_evaluate(yobs, event)
                                if mask is not None:
                                    xdata = xdata[ak.any(mask, axis=1)]
                                    ydata = ydata[ak.any(mask, axis=1)]
                                if ak.any(xdata) and ak.any(ydata):
                                    hist.from_array(xdata, ydata)
                                if hist.reach_limit:
                                    SSRLHistMaker.early_termination_counter += 1
                            elif hist.hist_type == "avg-graph":
                                xobs, yobs = hist.observable
                                try:
                                    xdata = ne_evaluate(xobs, event)
                                    ydata = ne_evaluate(yobs, event)
                                except IndexError:
                                    continue
                                if mask is not None:
                                    if mask.ndim == 1:
                                        xdata = xdata[mask]
                                        ydata = ydata[mask]
                                    else:
                                        xdata = xdata[ak.any(mask, axis=1)]
                                        ydata = ydata[ak.any(mask, axis=1)]
                                if ak.any(xdata) and ak.any(ydata):
                                    hist.from_array(xdata, ydata)
                            elif isinstance(hist, SSRLHisto1D):
                                hist.from_array(event, mask, w)
                            else:
                                obs = hist.observable[0]
                                data = ne_evaluate(obs, event)
                                if data.ndim != 1:
                                    data = ak.flatten(data)
                                if m_mask is not None:
                                    data = data[m_mask]
                                if ak.any(data):
                                    hist.from_array(data, w)

                    pbar_events.update(nevent)
                    if SSRLHistMaker.early_termination_counter > 20:
                        return p

        return p
