# -------------------------------------------------------------------
# extras interface
import MetXBase.extras
extras(m::AbstractHitOrDropSampler)::Dict = m.extras

# TODO: Revise this
# -------------------------------------------------------------------
# config interface
@extras_dict_interface AbstractHitOrDropSampler config

# -------------------------------------------------------------------
# state interface
@extras_dict_interface AbstractHitOrDropSampler state

# -------------------------------------------------------------------
# lep data
# NOTE: Do not interface with the net the data that will be in the constraints

import MetXBase.lepmodel
lepmodel(m::AbstractHitOrDropSampler) = m.lep