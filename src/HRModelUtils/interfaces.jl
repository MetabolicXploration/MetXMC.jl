# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::HRModel)::Dict = m.extras

# -------------------------------------------------------------------
# config interface
@extras_dict_interface HRModel config

# -------------------------------------------------------------------
# state interface
@extras_dict_interface HRModel state

# -------------------------------------------------------------------
# net interface
# NOTE: Do not interface with the net the data that will be in the constraints

# net data
@extras_val_interface HRModel metnet MetNet