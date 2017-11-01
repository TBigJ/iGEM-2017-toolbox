import numpy as np
from literature_parameters import literature_parameter_set
P_DEFAULT = literature_parameter_set(1)

def ode_function_actual_system(X, p = P_DEFAULT) :
    degradation_dynamics_n  = p ['P_tot'] * (p['k'] + p['delta']) * X / p['K']
    degradation_dynamics_d  = 1 + np.sum( X / p['K'])
    retroactivity_n         = (\
                                   1 +
                                   np.sum (X / p['K'])\
                                   - X / p['K']\
                              ) *\
                              p['P_tot'] / p['K']   
    retroactivity_d         = (1+np.sum ( X / p['K']))**2
    dXdt    =   ( p['H'] - degradation_dynamics_n / degradation_dynamics_d - p['delta'] * X ) /\
                (1+ retroactivity_n / retroactivity_d)
    return dXdt

def ode_function_no_retro(X, p = P_DEFAULT) :
    degradation_dynamics_n  = p ['P_tot'] * (p['k'] + p['delta']) * X / p['K']
    degradation_dynamics_d  = 1 + np.sum( X / p['K'])
    dXdt    =   ( p['H'] - degradation_dynamics_n / degradation_dynamics_d - p['delta'] * X )
    return dXdt

def ode_function_dilution_only(X, p = P_DEFAULT):
    return p['H']-X*p['delta']

def ode_function_no_saturation(X, p = P_DEFAULT):
    production  = p['H']
    proportionate_degradation = p['P_tot']*(p['k']+p['delta'])
    dilution    = p['delta']
    return production - proportionate_degradation * X - dilution * X

def ode_function_full_saturation(X, p = P_DEFAULT):
    production  = p['H']
    constant_degradation = p['P_tot']*(p['k']+p['delta'])/p['K']
    dilution = p['delta']
    return production - constant_degradation - dilution * X

ode_dict = {
           'actual':ode_function_actual_system,
           'no_retro':ode_function_no_retro,
           'dilution_only':ode_function_dilution_only,
           'no_saturation':ode_function_no_saturation,
            'full_saturation':ode_function_full_saturation
           }

def get_ode(system_type, parameters):
    return lambda t, Y: ode_dict[system_type](Y, parameters)