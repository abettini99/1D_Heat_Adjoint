import struct
import numpy as np
import numpy.typing as npt

## @brief Reads output of main.cpp executable, a binary data file, and outputs 2d arrays.
#
#  @details
#  The binary file holds data in the form:
#
#  u32  u32    f32  f32  f32    f32  f32  f32    ...   f32  f32  f32
#  jmax imax   x[0] y[0] u[0]   x[1] y[1] u[1]   ...   x[K] y[K] u[K],
#
#  where K = (jmax-1)*(imax-1). Reshapes x,y,u to (jmax,imax)
#
#  @param fileName Name of the binary file to read
#
#  @return x 2D numpy array of data x-positions
#  @return y 2D numpy array of data y-positions
#  @return u 2D numpy array of data values
def read(fileName: str) -> tuple[npt.NDArray[np.float32], 
                                 npt.NDArray[np.float32],
                                 npt.NDArray[np.float32]]:

    # Open file
    with open(fileName, mode='rb') as file: # b is important -> binary
        fileContent = file.read()

    ## ========= ##
    ## Read Data ##
    ## ========= ##
    # see link for more details on data types: https://docs.python.org/3/library/struct.html
    # Read jmax, imax as u32
    jmax = struct.unpack("I", fileContent[0:4])[0]
    imax = struct.unpack("I", fileContent[4:8])[0]

    # Read rest of data as f32, reshape it to the correct size.
    data = list( struct.unpack('f'*((len(fileContent) - 8) // 4), fileContent[8:]) )
    data = np.asarray(data, dtype=np.float32)
    data = data.reshape((-1,3))

    return  data[:,0].reshape((jmax,imax)), \
            data[:,1].reshape((jmax,imax)), \
            data[:,2].reshape((jmax,imax))  # (x,y,u)
