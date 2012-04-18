'''
Created on Mar 22, 2012

@author: ckernan
'''

import numpy as np

class MatrixWrapper(object):
    
    def __init__(self, columns, initial=100):
        """
        If initial is an integer, then it will create an empty matrix with 
        initial number of rows. If it is a numpy array, then it will just wrap
        that array.
        """
        self.columns = columns
        self._columnsDict = dict((col, index) 
                                for index, col in enumerate(columns))
        self._numCols = len(self.columns)
        
        if isinstance(initial, int):
            self._numRows = initial
            self._count = 0
            self.array = np.empty((initial, self._numCols), dtype='float')
        else:
            rows, cols = initial.shape
            assert cols == self._numCols
            
            self._numRows = rows
            self._count = rows
            self.array = initial
        
    @property
    def shape(self):
        return self._numRows, self._numCols
        
    def append(self, *data):
        """ Data needs to be in the same order as columns. """
        
        self.array[self._count, :] = data
        self._count += 1
        
        if self._count >= self._numRows:
            self._numRows *= 2
            self.array.resize(self.shape)
            
    def compact(self):
        self._numRows = self._count
        self.array.resize(self.shape)
        
    def selectColumns(self, *columns):
        newArray = np.empty((self._numRows, len(columns)), dtype='float')
        for i, column in enumerate(columns):
            newArray[:, i] = self.array[:, self._columnsDict[column]]
        return newArray
    
    def merge(self, matrix):
        assert self.columns == matrix.columns
        # TODO this erases subclass info
        
        newArray = np.append(self.array, matrix.array, 0)
        return MatrixWrapper(self.columns, newArray)
    
    def rows(self):
        for row in range(self._numRows):
            yield self.array[row, :]
    
    def __getitem__(self, index):
        rows, columns = index
        
        if isinstance(columns, tuple):
            ndmin = 1 if isinstance(rows, int) else 2
            first = True
            for column in columns:
                colIndex = self._columnsDict[column]
                data = np.array(self.array[rows, colIndex], ndmin=ndmin)
                
                if first:
                    items = data
                    first = False
                else:
                    items = np.append(items, data, 0)
            
            return items.transpose()
        elif columns == slice(None):
            return self.array[rows, :]
        else:
            return self.array[rows, self._columnsDict[columns]]
        
if __name__ == '__main__':
    print "YAY"
    
    class Test(MatrixWrapper):
        
        def __init__(self):
            a = np.arange(9).reshape((3, 3))
            print a
            MatrixWrapper.__init__(self, ["x", "y", "z"], a)
    
    t = Test()
    print t[:, "x"]
    
