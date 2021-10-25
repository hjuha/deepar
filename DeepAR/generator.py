import random

group, n = input().split()[:3]
n = int(n)
random.seed(n)

matrix = [[0 for _ in range(n)] for _ in range(n)]

print(n)

if group == "uniform":
	for y in range(n):
		for x in range(n):
			matrix[y][x] = random.random()
		print(" ".join([str(z) for z in matrix[y]]))
elif group == "block_diagonal":
	k = 5
	for y in range(n):
		for x in range(n):
			if (y // k == x // k):
				matrix[y][x] = random.random()
		print(" ".join([str(z) for z in matrix[y]]))
elif group == "staircase":
            for y in range(n):
                    for x in range(n):
                            if (x <= y + 1):
                                    matrix[y][x] = 1
                            else:
                                    matrix[y][x] = 0
                    print(" ".join([str(z) for z in matrix[y]]))