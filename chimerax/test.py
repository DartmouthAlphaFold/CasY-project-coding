import numpy as np
import pandas as pd
import time
from functools import wraps
from pvclust import PvClust

def timer(progress_func):
    '''
    Measure how much time each function executed
    '''
    @wraps(progress_func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = progress_func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Step {progress_func.__name__} took {execution_time:.3f} seconds")
        return result
    return wrapper

@timer
def main():
    data_url = "http://lib.stat.cmu.edu/datasets/boston"
    raw_df = pd.read_csv(data_url, sep="\s+", skiprows=22, header=None)
    X = np.hstack([raw_df.values[::2, :], raw_df.values[1::2, :2]])
    X = pd.DataFrame(X)
    pv = PvClust(X, method="ward", metric="euclidean", nboot=100, parallel=True)
    pv.plot(param_display='AU', sig_level = 95, orientation = "top")
    pv.print_result(digits = 5)
    pv.seplot(annotate=True)
    
    
if __name__ == "__main__":
    main()