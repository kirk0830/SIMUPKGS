# this function provides periodic boundary condition relationg operation
import numpy as np

def cell_operator(a, b, c, alpha, beta, gamma):

    """
    This function will return an numpy incidence in 2 dimensions, i.e. a matrix operator\n
    alpha: angle between a and b\n
    beta: angle between b and c\n
    gamma: angle between c and a\n
    """
    a = float(a)
    b = float(b)
    c = float(c)
    alpha_r = float(alpha/180) * np.pi
    beta_r = float(beta/180) * np.pi
    gamma_r = float(gamma/180) * np.pi
    a1 = np.round(a, 10)
    a_vec = [a1, 0, 0]
    b1 = np.round(b*np.cos(alpha_r), 10)
    b2 = np.round(b*np.sin(alpha_r), 10)
    b_vec = [b1, b2, 0]
    c1 = np.round(c*np.cos(gamma_r), 10)
    c2 = np.round(c*(np.cos(beta_r)-np.cos(alpha_r)*np.cos(gamma_r))/np.sin(alpha_r), 10)
    c3 = np.round(np.sqrt(c**2 - c1**2 - c2**2), 10)
    c_vec = [c1, c2, c3]
    # warning: for matrix calculation, all matrix should be standardized! 
    O = np.reshape([a_vec, b_vec, c_vec], (3,3)).T
    return O

# coordinates conversion:
# matrix elements in real space:
# a_vec_1 a_vec_2 a_vec_3        i_crys        x
# b_vec_1 b_vec_2 b_vec_3   *    j_crys    =   y
# c_vec_1 c_vec_2 c_vec_3        k_crys        z
# in succint form:
# O*r_crys = r_Cart

def convert2cart(O_in_cart, vector_crys):

    """
    Warning: O_in_cart must have standard numpy ndarray data type.
    """
    vector_crys = np.reshape(vector_crys, (3,1))
    return np.matmul(O_in_cart, vector_crys)

def convert2crys(O_in_cart, vector_cart):

    """
    Warning: O_in_cart must have standard numpy ndarray data type.
    """
    vector_cart = np.reshape(vector_cart, (3,1))
    return np.matmul(np.linalg.inv(O_in_cart), vector_cart)

def vectorPbcCorr(O_in_cart, displacement):
    
    # methodology: if projection of displacement is larger than half of a crystal vector,
    #              then substract it to get the new displacement vector
    # Also plz note that O_in_cart must be numpy.ndarray type data
    crys_a = O_in_cart[:][0]
    crys_b = O_in_cart[:][1]
    crys_c = O_in_cart[:][2]
    norm_a = np.linalg.norm(crys_a)
    norm_b = np.linalg.norm(crys_b)
    norm_c = np.linalg.norm(crys_c)
    proj_a = np.dot(crys_a, displacement)/norm_a
    proj_b = np.dot(crys_b, displacement)/norm_b
    proj_c = np.dot(crys_c, displacement)/norm_c
    if proj_a > norm_a/2:
        for iaxis in range(3):
            displacement[iaxis] -= crys_a[iaxis]
    if proj_b > norm_b/2:
        for iaxis in range(3):
            displacement[iaxis] -= crys_b[iaxis]
    if proj_c > norm_c/2:
        for iaxis in range(3):
            displacement[iaxis] -= crys_c[iaxis]

    return displacement