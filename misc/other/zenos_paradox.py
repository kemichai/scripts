"""
...

Date: Mar 2023
Author: KM
"""

print("Suppose someone wishes to step to a door...")
x = float(input("Enter the distance to that door (let's say in meters): "))

cov_dist = 0
rem_dist = x
i = 1
while cov_dist < x:
# while i < 5:
    step = rem_dist/2
    cov_dist += step
    rem_dist = x - cov_dist
    # print(f"Covered distance: {cov_dist}, Remaining distance: {rem_dist}")
    ratio = cov_dist / x
    rat = ratio.as_integer_ratio()
    print(f"After step {i} {rat[0]}/{rat[1]} distance covered...")
    i += 1


