import bempp.api 
import numpy as np
import time
import scipy.linalg

def rescale(A, d1, d2):
    """Rescale the 2x2 block operator matrix A"""
    
    A[0, 1] = A[0, 1] * (d2 / d1)
    A[1, 0] = A[1, 0] * (d1 / d2)
    
    return A

def rotate_x(x,y,z,theta):
    c = np.cos(theta)
    s = np.sin(theta)
    R_x = np.array([[1,0,0],[0,c, -s], [0,s, c]])
    v = np.array([x,y,z])
    rotated_v = R_x.dot(v)
    return rotated_v

def rotate_y(x,y,z,theta):
    c = np.cos(theta)
    s = np.sin(theta)
    R_y = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    v = np.array([x,y,z])
    rotated_v = R_y.dot(v)
    return rotated_v

def rotate_z(x,y,z,theta):
    c = np.cos(theta)
    s = np.sin(theta)
    R_z = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    v = np.array([x,y,z])
    rotated_v = R_z.dot(v)
    return rotated_v

def translation(x,y,z,d):
    x += d[0]
    y += d[1]
    z += d[2]
    return (x,y,z)

def PMCHWT_operator(grids, k_ext, k_int, mu_ext, mu_int, parameters = None):
    number_of_scatterers = len(grids)
    result = bempp.api.BlockedOperator(2*number_of_scatterers, 2*number_of_scatterers)
    
    rwg_space = [bempp.api.function_space(grid, 'RWG', 0) for grid in grids]
    snc_space = [bempp.api.function_space(grid, 'SNC', 0) for grid in grids]
    
    interior_operators = []
    exterior_operators = []
    identity_operators = []
    
    for i in range(number_of_scatterers):
        Ai = bempp.api.BlockedOperator(2,2)
        Ae = bempp.api.BlockedOperator(2,2)

        magnetic_ext = bempp.api.operators.boundary.maxwell.magnetic_field(domain=rwg_space[i], range_=rwg_space[i], 
                                                                           dual_to_range=snc_space[i], 
                                                                           wave_number=k_ext, parameters = parameters)
        electric_ext = bempp.api.operators.boundary.maxwell.electric_field(domain=rwg_space[i], range_=rwg_space[i],
                                                                          dual_to_range = snc_space[i], wave_number = k_ext,
                                                                          parameters = parameters)

        magnetic_int = bempp.api.operators.boundary.maxwell.magnetic_field(domain=rwg_space[i], range_=rwg_space[i], 
                                                                           dual_to_range=snc_space[i], wave_number=k_int[i],
                                                                          parameters = parameters)
        electric_int = bempp.api.operators.boundary.maxwell.electric_field(domain=rwg_space[i], range_=rwg_space[i],
                                                                          dual_to_range = snc_space[i], wave_number = k_int[i],
                                                                           parameters=parameters)

        op_ident = bempp.api.assembly.BlockedOperator(2,2)
        op_identity = bempp.api.operators.boundary.sparse.identity(rwg_space[i],rwg_space[i],snc_space[i])
        op_ident[0,0] = op_identity
        op_ident[1,1] = op_identity
        identity_operators.append(op_ident)

        Ai[0,0] = magnetic_int
        Ai[0,1] = electric_int
        Ai[1,0] = -1* electric_int
        Ai[1,1] = magnetic_int

        Ai = rescale(Ai, k_int[i], mu_int[i])
        interior_operators.append(Ai)

        Ae[0,0] = magnetic_ext
        Ae[0,1] = electric_ext
        Ae[1,0] = -1* electric_ext
        Ae[1,1] = magnetic_ext

        Ae = rescale(Ae, k_ext, mu_ext)
        exterior_operators.append(Ae)

    filter_operators = []
    transfer_operators = np.empty((number_of_scatterers, number_of_scatterers), dtype=np.object)
    for i in range(number_of_scatterers):
        filter_operators.append(0.5*identity_operators[i] - interior_operators[i])
        for j in range(number_of_scatterers):
            if i==j:
                element = interior_operators[i]+exterior_operators[j]
            else:
                Aij = bempp.api.BlockedOperator(2,2)

                magnetic_ij = bempp.api.operators.boundary.maxwell.magnetic_field(domain = rwg_space[j], range_ = rwg_space[i],
                                                                                 dual_to_range = snc_space[i], wave_number = k_ext,
                                                                                 parameters = parameters)
                electric_ij = bempp.api.operators.boundary.maxwell.electric_field(domain = rwg_space[j], range_ = rwg_space[i],
                                                                                 dual_to_range = snc_space[i], wave_number = k_ext,
                                                                                 parameters = parameters)

                Aij[0,0] = magnetic_ij
                Aij[0,1] = electric_ij
                Aij[1,0] = -1* electric_ij
                Aij[1,1] = magnetic_ij

                Aij = rescale(Aij, k_ext, mu_ext)

                transfer_operators[i,j] = Aij
                element = transfer_operators[i,j]

            result[2 * i, 2 * j] = element[0, 0]
            result[2 * i, 2 * j + 1] = element[0, 1]
            result[2 * i + 1, 2 * j] = element[1, 0]
            result[2 * i + 1, 2 * j + 1] = element[1, 1] 
            
    return [result, filter_operators]

def PMCHWT_preconditioner(grids, k_ext, k_int, mu_ext, mu_int, parameters = None, full = False, electric_only = False):
    number_of_scatterers = len(grids)
    result = bempp.api.BlockedOperator(2*number_of_scatterers, 2*number_of_scatterers)

    bc_space = [bempp.api.function_space(grid, 'BC', 0) for grid in grids]
    rbc_space = [bempp.api.function_space(grid, 'RBC', 0) for grid in grids]
    b_rwg_space = [bempp.api.function_space(grid, "B-RWG", 0) for grid in grids]
    
          
    for i in range(number_of_scatterers):
        Ai_pre = bempp.api.BlockedOperator(2,2)
        Ae_pre = bempp.api.BlockedOperator(2,2)

        magnetic_ext = bempp.api.operators.boundary.maxwell.magnetic_field(domain = bc_space[i],range_= b_rwg_space[i], 
                                                                       dual_to_range = rbc_space[i], wave_number = k_ext,
                                                                      parameters = parameters)
        electric_ext = bempp.api.operators.boundary.maxwell.electric_field(domain = bc_space[i],range_= b_rwg_space[i], 
                                                                       dual_to_range = rbc_space[i], wave_number = k_ext,
                                                                      parameters = parameters)

        magnetic_int = bempp.api.operators.boundary.maxwell.magnetic_field(domain = bc_space[i],range_= b_rwg_space[i], 
                                                                       dual_to_range = rbc_space[i], wave_number = k_int[i],
                                                                      parameters = parameters)
        electric_int = bempp.api.operators.boundary.maxwell.electric_field(domain = bc_space[i],range_= b_rwg_space[i], 
                                                                       dual_to_range = rbc_space[i], wave_number = k_int[i],
                                                                      parameters = parameters)

        Ai_pre[0,0] = magnetic_int
        Ai_pre[0,1] = electric_int
        Ai_pre[1,0] = -1* electric_int
        Ai_pre[1,1] = magnetic_int

        Ai_pre = rescale(Ai_pre, k_int[i], mu_int[i])

        Ae_pre[0,0] = magnetic_ext
        Ae_pre[0,1] = electric_ext
        Ae_pre[1,0] = -1* electric_ext
        Ae_pre[1,1] = magnetic_ext

        Ae_pre = rescale(Ae_pre, k_ext, mu_ext)

        
        if full == True:
            element = Ai_pre + Ae_pre
        else:
            element = Ai_pre

        if electric_only == True:
            result[2*i, 2*i+1] = element[0,1]
            result[2*i+1, 2*i] = element[1,0]
        else:
            result[2 * i, 2 * i] = element[0, 0]
            result[2 * i, 2 * i + 1] = element[0, 1]
            result[2 * i + 1, 2 * i] = element[1, 0]
            result[2 * i + 1, 2 * i + 1] = element[1, 1] 
        
    return result


def mass_matrix_BC_SNC(grids):
    number_of_scatterers = len(grids)
    
    bc_space = [bempp.api.function_space(grid, 'BC', 0) for grid in grids]

    b_rwg_space = [bempp.api.function_space(grid, "B-RWG", 0) for grid in grids]
    b_snc_space = [bempp.api.function_space(grid, "B-SNC", 0) for grid in grids]
    
    result = np.empty((2*number_of_scatterers,2*number_of_scatterers), dtype = 'O')

    for i in range(number_of_scatterers):
        ident = bempp.api.operators.boundary.sparse.identity(bc_space[i], b_rwg_space[i], b_snc_space[i])
        inv_ident = bempp.api.assembly.InverseSparseDiscreteBoundaryOperator(ident.weak_form())
        result[2*i,2*i] = inv_ident
        result[2*i+1,2*i+1] = inv_ident

    result = bempp.api.assembly.BlockedDiscreteOperator(result)
    return result





