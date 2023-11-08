

# Our model is composed of the following parts:
#
# Forces 
#   internal forces acting on a single agent
#   local interaction forces acting between closeby agents
#   graph interaction forces acting on connected agents  
#
# Events
#   changes of one agent's parameters 
#   changes affecting the whole system
#
# Enviroment
#   the enviroment might use a model itself, which is assumed to be 
#   updated between timesteps (e.g. reaction-diffusion model for signals)
#   and interaction with agents is resolved once per time-step 



A = Agent()
B = Agent()