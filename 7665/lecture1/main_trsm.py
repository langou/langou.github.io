import time
import numpy as np
import copy

def forward__substitution__scalar( A, b ):
  n = b.shape[0]
  p = b.shape[1]
  for i in range(0,n):
    for j in range(0,i):
      for k in range(0,p):
        b[i,k] = b[i,k] - A[i,j] * b[j,k]
    for k in range(0,p):
      b[i,k] = b[i,k] / A[i,i]

def forward__substitution__vector( A, b ):
  n = b.shape[0]
  k = b.shape[1]
  for i in range(0,n):
    for j in range(0,i):
      b[i,0:k] = b[i,0:k] - A[i,j] * b[j,0:k]
    b[i,0:k] = b[i,0:k] / A[i,i]

def forward__substitution__matrix(A: np.ndarray, b: np.ndarray) -> np.ndarray:
  n = b.shape[0]
  k = b.shape[1]
  for i in range(n):
      b[i,0:k] = (b[i,0:k] - A[i, :i] @ b[:i,0:k]) / L[i, i]

#############################

n = 4000
k = 200

L = np.tril(np.random.randn(n, n), -1)*1e-6
np.fill_diagonal(L, 1)

b = np.random.randn(n,k)

#############################

x = np.copy(b)
elapsed_time =- time.time()
forward__substitution__matrix(L, x)
elapsed_time += time.time()

residual = np.linalg.norm(L @ x - b) / np.linalg.norm(b) 
print(f"Matrix size: {n}x{n}")
print(f"Execution time: {elapsed_time:.6f} seconds")
print(f"Execution time: {n*n*k/elapsed_time*1e-9:.3f} GFlop/sec")
print(f"Residual ||Lx - b||: {residual:.2e}\n")

#############################

x = np.copy(b)

elapsed_time =- time.time()
forward__substitution__vector(L, x)
elapsed_time += time.time()

residual = np.linalg.norm(L @ x - b)
print(f"Matrix size: {n}x{n}")
print(f"Execution time: {elapsed_time:.6f} seconds")
print(f"Execution time: {n*n*k/elapsed_time*1e-9:.3f} GFlop/sec")
print(f"Residual ||Lx - b||: {residual:.2e}\n")

#############################

x = np.copy(b)

elapsed_time =- time.time()
forward__substitution__scalar(L, x)
elapsed_time += time.time()

residual = np.linalg.norm(L @ x - b)
print(f"Matrix size: {n}x{n}")
print(f"Execution time: {elapsed_time:.6f} seconds")
print(f"Execution time: {n*n*k/elapsed_time*1e-9:.3f} GFlop/sec")
print(f"Residual ||Lx - b||: {residual:.2e}\n")

#############################

x = np.copy(b)

elapsed_time =- time.time()
x = np.linalg.solve(L, b)
elapsed_time += time.time()

residual = np.linalg.norm(L @ x - b)
print(f"Matrix size: {n}x{n}")
print(f"Execution time: {elapsed_time:.6f} seconds")
print(f"Execution time: {n*n*k/elapsed_time*1e-9:.3f} GFlop/sec")
print(f"Residual ||Lx - b||: {residual:.2e}\n")

#############################

