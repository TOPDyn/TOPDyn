import os
import numpy as np

class SaveData():
    def __init__(self, constr_func, max_iter, multi_obj_bool, save):

        self.constr_func = constr_func

        if save:
            folder_name = 'data'
            self.directory = os.path.join(os.path.dirname(__file__), folder_name)
            os.makedirs(self.directory, exist_ok=True)

            self.xnew_dir = os.path.join(self.directory, 'xnew')
            os.makedirs(self.xnew_dir, exist_ok=True)

            self.save_fval = np.empty((max_iter + 1, len(self.constr_func)))
            self.save_f0val1 = np.empty(max_iter + 1)
            self.save_fvirg1 = np.empty(max_iter + 1)

            if multi_obj_bool:
                self.save_f0val2 = np.empty(max_iter + 1)
                self.save_fvirg2 = np.empty(max_iter + 1)

    def update_save_f0val1(self, outit, f0val1, fvirg1):
        self.save_fvirg1[outit] = fvirg1
        self.save_f0val1[outit] = f0val1

    def update_save_f0val2(self, outit, f0val2, fvirg2):
        self.save_fvirg2[outit] = fvirg2
        self.save_f0val2[outit] = f0val2

    def update_save_fval(self, outit, fval):
        self.save_fval[outit, :] = fval[:, 0]    

    def save_xval(self, outit, xval):
        np.savetxt(os.path.join(self.xnew_dir, 'xnew'+'_'+str(outit)+'.txt'), xval)

    def _create_header(self, multiobj_bool):
        """ Creates header to save data.

        Args:
            multiobj_bool (:obj:`bool`): True if multiobjective function is used.
        
        Returns:
            Header.
        """
        if multiobj_bool:
            header = "iter,f0val,fvirg,f0val2,fvirg2"
        else:
            header = "iter,f0val,fvirg"

        for func in self.constr_func:
            header += ',' + func
        return header

    def save_data(self, multiobj_bool, outit, f_original=None, f_optimized=None):
        if multiobj_bool:
            ind_fval = 5
            aux_multi = 2
        else:
            ind_fval = 3
            aux_multi = 0
         
        data = np.empty((outit + 1, 3 + len(self.constr_func) + aux_multi))
        data[:, 0] = np.arange(1, outit+2)

        data[:, 1] = self.save_f0val1[:outit + 1]
        data[:, 2] = self.save_fvirg1[:outit + 1]

        if multiobj_bool:
            data[:, 3] = self.save_f0val2[:outit + 1]
            data[:, 4] = self.save_fvirg2[:outit + 1]

        data[:, ind_fval:] = self.save_fval[:outit + 1,:]

        header = self._create_header(multiobj_bool)
        np.savetxt(os.path.join(self.directory, 'functions.txt'), data, delimiter=",", header=header, comments='')

        if f_original is not None:
            np.savetxt(os.path.join(self.directory, 'frequency_rsp.txt'), np.column_stack((f_original, f_optimized)), delimiter=",", header="orig,opt", comments='')