# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::AbstractHitOrDropSampler)::Dict = m.extras

# -------------------------------------------------------------------
# config interface
@extras_dict_interface AbstractHitOrDropSampler config

# -------------------------------------------------------------------
# state interface
@extras_dict_interface AbstractHitOrDropSampler state

# -------------------------------------------------------------------
# net interface
# NOTE: Do not interface with the net the data that will be in the constraints

# net data
@extras_val_interface AbstractHitOrDropSampler metnet MetNet