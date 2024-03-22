import lazy_loader as lazy

__getattr__, __dir__, __all__ = lazy.attach(
    __name__,
    # submodules=['strategies', 'histManipulate', 'run_HistMaker', 'run_PlotMaker'],
    submod_attrs={
        'version': ['__version__'],
        'histmaker': ['Graph', 'SSRLHisto1D', 'AvgGraph', 'SSRLHistMaker'],
    },
)
