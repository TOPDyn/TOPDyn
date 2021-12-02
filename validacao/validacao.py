import numpy as np
import os

def validate_xnew(path_correto, path_check):
    files_correto = os.listdir(path_correto)
    for file in files_correto:
        xnew_correto = np.round(np.loadtxt(os.path.join(path_correto, file)), 6)
        xnew_check = np.round(np.loadtxt(os.path.join(path_check, file)), 6)
        if not np.array_equal(xnew_correto, xnew_check):
            print("Damn it! - validate_xnew")
            print(file+"\n")
            break

def validade_funcs(path_correto, path_check):
    df_correto = np.round(np.loadtxt(path_correto, delimiter=',', skiprows=1),7)
    df_check = np.round(np.loadtxt(path_check, delimiter=',', skiprows=1),7)
    if not np.array_equal(df_correto, df_check):
        print("Damn it! - validade_funcs")
        print(path_check+"\n")

def validade_freqrsp(path_correto, path_check):
    freq_correto = np.loadtxt(path_correto, delimiter=',', skiprows=1, dtype=np.complex)
    freq_check = np.loadtxt(path_check, delimiter=',', skiprows=1, dtype=np.complex)
    if not np.array_equal(freq_correto, freq_check):
        print("Damn it! - validade_freqrsp")
        print(path_check+"\n")

correto_xnew = "D:\\TOPDyn\\validacao\\data_correto_gcmma"

check_xnew = "D:\\TOPDyn\\data"

validate_xnew(os.path.join(correto_xnew, "xnew"), os.path.join(check_xnew, "xnew"))

validade_funcs(os.path.join(correto_xnew, "functions.txt"), os.path.join(check_xnew, "functions.txt"))

validade_freqrsp(os.path.join(correto_xnew, "frequency_rsp.txt"), os.path.join(check_xnew, "frequency_rsp.txt"))


