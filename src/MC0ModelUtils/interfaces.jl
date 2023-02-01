# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::MC0Model)::Dict = m.extras

# -------------------------------------------------------------------
# config interface
@extras_dict_interface MC0Model config

# -------------------------------------------------------------------
# state interface
@extras_dict_interface MC0Model state

# -------------------------------------------------------------------
# net interface
# NOTE: Do not interface with the net the data that will be in the constraints

# net data
@extras_val_interface MC0Model metnet MetNet