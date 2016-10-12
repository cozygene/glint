import logging
import cvxopt
from utils import common
from numpy import loadtxt, zeros, matrix, array, savetxt, column_stack, vstack, double
from module import  Module



class Houseman(Module):
    
    def __init__(self, methylation_data, reference_file, outputfile = None):
        """
         reference_file is path to a .txt file of sites by cell type refernece matrix, including row and column headers.
        """
        self.meth_data = methylation_data

        # Load and extract the reference - cell-type specific methylation levels of a group of reference sites
        logging.info("Loading houseman refernece file %s..." % reference_file.name)
        ref_data, self.names, ref_ids = common.load_data_file(reference_file.name, 2)
        if ref_data is None or ref_ids is None:
            common.terminate("There is a problem with the format of the reference file '%s'." % reference_file.name)

        self.components = self.houseman(ref_data, ref_ids)
        self.save(outputfile)
        self.outputfile = outputfile

    def houseman(self, ref_data, ref_ids):
        """
         function for estimating cell counts
        """
        cvxopt.solvers.options['show_progress'] = False

        # keep only site that are in the reference list
        self.meth_data.include(ref_ids) 

        # use only the reference sites that are available in the data, and sort the reference data in the same order
        d = {}
        counter = 0
        for i in ref_ids:
            d[i] = counter
            counter += 1

        keep = []
        for i in self.meth_data.cpgnames:
            keep.append(d[i])

        ref_ids_final = ref_ids[keep]
        ref_data_final = ref_data[keep,:]

        n = len(self.meth_data.data[0])
        m_cvxopt = cvxopt.matrix(ref_data_final.astype(double)) #CVX only accepts doubles, not floats
        o_cvxopt = cvxopt.matrix(self.meth_data.data.astype(double))
        m, k = m_cvxopt.size
        R_est = zeros((n, k)) 
        for i in range(n): 
            Q =  cvxopt.matrix(2*m_cvxopt.T*m_cvxopt, tc ='d')
            p = (-2*o_cvxopt.T*m_cvxopt).T[:,i]
            I = cvxopt.matrix(0.0, (k,k))
            I[::k+1] = -1.0
            G = I
            h = cvxopt.matrix(k*[0.0])
            b = cvxopt.matrix(1*[1.0])
            A = cvxopt.matrix(k*[1.0]).T
            sol= cvxopt.solvers.qp(Q, p, G, h, A, b)
            r_hat = sol['x']
            R_est[i,:] = matrix(array(r_hat)).T
        
        return R_est

    def save(self, outputfile):
        if outputfile:
            logging.info("Saving Houseman components to %s..." % outputfile)
            components_output = column_stack((self.meth_data.samples_ids, self.components))
            header = ["ID"] + list(self.names)
            components_output = vstack((header, components_output))
            savetxt(outputfile+"S", components_output, fmt='%s') #output looks betther with space delimiter and not tab
            


