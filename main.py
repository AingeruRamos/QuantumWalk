import sys
import time
import math
import cmath

import numpy as np

from scipy.optimize import minimize
from scipy.linalg import fractional_matrix_power

#==================================================================#
# SHIFT MATRIX FUNCTIONS                                           #
#==================================================================#

#
# Gets the total number of qbits requested by variables
#
def get_vars_q_size(var_list):
  return np.sum([*map(lambda x: int(''.join(x[1:])), var_list)])

#
# Calculates the total operations needed for the variables
#
def get_n_ops(var_list):
  def l(var_info):
    t, *q_size = var_info
    q_size = int(''.join(q_size))
    return 2 if t == 'c' else pow(2, q_size)-1

  return np.sum([*map(l, var_list)])

#
# Returns the list of operations of a continuous variable
#
def c_ops(q_size):
  ident = np.identity(pow(2, q_size))
  return [np.roll(ident, -1, axis=0),
          np.roll(ident, 1, axis=0)]

#
# Returns the list of operations of a discrete variable
#
def d_ops(q_size):
  statespace_size = pow(2, q_size)
  ident = np.identity(statespace_size)
  return [np.roll(ident, i, axis=0) for i in range(1, statespace_size)]

#
# Return the list of operations of any type of variable.
# Acts like a wrapper for ''c_ops'' and ''d_ops'' functions.
#
def get_var_subops(var_info):
  t, *q_size = var_info
  q_size = int(''.join(q_size))

  if t == 'c': return c_ops(q_size)
  elif t == 'd': return d_ops(q_size)

#==================================================================#

#
# Given a list of variables, returns the matrix S representing the all posible movements.
#
def compose_operation(var_list):

  pos_statespace_size = pow(2, get_vars_q_size(var_list))

  n_ops = get_n_ops(var_list)
  dice_statespace_size = pow(2, math.ceil(math.log2(n_ops)))

  circuit_statespace_size = pos_statespace_size*dice_statespace_size

  final_op = np.identity(circuit_statespace_size)

  corner_indexes = [0, pos_statespace_size]

  for i, var_info in enumerate(var_list):
    for subop in get_var_subops(var_info):
      composed_subop = [1]
      for j, var_info in enumerate(var_list):
        q_size = int(''.join(var_info[1:]))
        if i == j: composed_subop = np.kron(composed_subop, subop)
        elif i != j: composed_subop = np.kron(composed_subop, np.identity(pow(2, q_size)))

      final_op[corner_indexes[0]:corner_indexes[1], corner_indexes[0]:corner_indexes[1]] = composed_subop
      corner_indexes[0] += pos_statespace_size
      corner_indexes[1] += pos_statespace_size

  return final_op

#==================================================================#
# UNITARY MATRIX EQUATION                                          #
#==================================================================#

#
# Given 2 real vectors, builds a complex type vector
#
def get_complex_vector(v_real, v_imag):
  assert len(v_real) == len(v_imag)
  return np.array([complex(v_real[i], v_imag[i]) for i in range(len(v_real))])

#
# Calculates scalar*(v1*conjugate(v2))-b
# Used to compose matrix constraints
#
def get_2Vector_constraint(v1, v2, b, scalar=1):
  assert len(v1) == len(v2)

  c_val = 0+0j
  for i in range(len(v1)):
    c_val += scalar*v1[i]*v2[i].conjugate()
  
  # Bastante fumada
  return (abs(c_val.real) - b) + abs(c_val.imag)

#==================================================================#

#
# Compose the constraint to get a unitary matrix
#
def get_unitary_constraints(m_size):
  constraints = []

  def constraint_f_factory(i, j):
    def constraint_f(x):
        v1_complex = get_complex_vector(x[i*m_size:i*m_size+m_size], x[2*i*m_size:2*i*m_size+m_size])
        v2_complex = get_complex_vector(x[j*m_size:j*m_size+m_size], x[2*j*m_size:2*j*m_size+m_size])
        return get_2Vector_constraint(v1_complex, v2_complex, int(i==j))
    return constraint_f

  for i in range(m_size):
    for j in range(m_size):
      constraints.append(
        {'type': 'ineq', 'fun': constraint_f_factory(i, j)}
      )
  return constraints

#
# Compose the objective function to get the operation we want
#
def optimize_func(x, V_end, m_size):
  M = get_complex_vector(x[0:m_size*m_size], x[m_size*m_size:])
  M = M.reshape((m_size, m_size))

  V_start = m_size * [math.sqrt(1.0/(m_size))]
  V_start = np.array(V_start)

  V_process = V_end-np.matmul(M, V_start)
  return np.linalg.norm(V_process)

#==================================================================#
# MAIN                                                             #
#==================================================================#

# Initialize constants
#==================================================================#
C_Q_SIZE = int(sys.argv[1])
D_Q_SIZE = int(sys.argv[2])
R_INDEX = int(sys.argv[3])-1
FILENAME = sys.argv[4]

ITERATIONS = 1

VARS = [f'c{C_Q_SIZE}', f'd{D_Q_SIZE}']
if D_Q_SIZE == 0:
  del VARS[1]

# Get S Matrix
#==================================================================#

time_start = time.time()

S = compose_operation(VARS)
S_INV = np.transpose(S)

# Set optimization problem
#==================================================================#

circuit_statespace_size = S.shape[0]
position_statespace_size = pow(2, get_vars_q_size(VARS))
dice_statespace_size = int(circuit_statespace_size/position_statespace_size) 

M_initial = 2*np.random.random(2*circuit_statespace_size*circuit_statespace_size)-1

# JUST FOR DEBUG
#==================================================================#
M_real = M_initial[0:circuit_statespace_size*circuit_statespace_size]
M_imag = M_initial[circuit_statespace_size*circuit_statespace_size:]
M = get_complex_vector(M_real, M_imag)
M = M.reshape((circuit_statespace_size, circuit_statespace_size))
print(M)
print()
#==================================================================#

V_end = circuit_statespace_size * [math.sqrt(1.0/circuit_statespace_size)]
if R_INDEX > -1:
  V_end = position_statespace_size * [0.0]
  V_end[0] = m.sqrt(1.0/dice_statespace_size)
  V_end = np.roll(V_end, R_INDEX).tolist()
  V_end = dice_statespace_size * V_end
V_end = np.array(V_end)

# Execute optimizer
#==================================================================#

objective = lambda x: optimize_func(x, V_end, circuit_statespace_size)
x0 = M_initial
bnds = (2*circuit_statespace_size*circuit_statespace_size) * [(-1.0, 1.0)]
cons = get_unitary_constraints(circuit_statespace_size)

res = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints=cons, tol=1e-6)
print(res)

M_real = res.x[0:circuit_statespace_size*circuit_statespace_size]
M_imag = res.x[circuit_statespace_size*circuit_statespace_size:]
M = get_complex_vector(M_real, M_imag)
M = M.reshape((circuit_statespace_size, circuit_statespace_size))

# Get C matrix from M
#==================================================================#

M = fractional_matrix_power(M, 1/ITERATIONS)

C = np.matmul(S_INV, M)

time_end = time.time()