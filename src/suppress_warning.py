import warnings
from numba.core.errors import NumbaTypeSafetyWarning

warnings.filterwarnings('ignore', category=NumbaTypeSafetyWarning)
warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=FutureWarning)

# Add NumPy cast warning suppression
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in cast')