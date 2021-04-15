
# p = (1, 0, 0)
p = zeros(4)
p[2] = 1.0
p[3] = 0.0
p[4] = 0.0


# v = (0, 0, 1)
v = zeros(3)
v[3] = 1.0

# theta
# 右手系と逆になる
# Mz > 0 ならば負の方向に回転
theta = -pi/2

# q = (cos(t/2),  vx*sin(t/2),  vy*sin(t/2),  vz*sin(t/2))
q = zeros(4)

q[1] = cos(theta/2)
q[2] = v[1]*sin(theta/2)
q[3] = v[2]*sin(theta/2)
q[4] = v[3]*sin(theta/2)



A = zeros(3,3)

A[1,1] = 1 - 2*q[3]^2 - 2*q[4]^2
A[1,2] = 2*q[2]*q[3] + 2*q[1]*q[4]
A[1,3] = 2*q[2]*q[4] - 2*q[1]*q[3]

A[2,1] = 2*q[2]*q[3] - 2*q[1]*q[4]
A[2,2] = 1 - 2*q[4]^2 - 2*q[2]^2
A[2,3] = 2*q[3]*q[4] + 2*q[1]*q[2]

A[3,1] = 2*q[2]*q[4] + 2*q[1]*q[3]
A[3,2] = 2*q[3]*q[4] - 2*q[1]*q[2]
A[3,3] = 1 - 2*q[2]^2 - 2*q[3]^2


newx = A[1,1]*p[2] + A[1,2]*p[3] + A[1,3]*p[4]
newy = A[2,1]*p[2] + A[2,2]*p[3] + A[2,3]*p[4]
newz = A[3,1]*p[2] + A[3,2]*p[3] + A[3,3]*p[4]

println(newx)
println(newy)
println(newz)

p[2] = newx
p[3] = newy
p[4] = newz

newx = A[1,1]*p[2] + A[1,2]*p[3] + A[1,3]*p[4]
newy = A[2,1]*p[2] + A[2,2]*p[3] + A[2,3]*p[4]
newz = A[3,1]*p[2] + A[3,2]*p[3] + A[3,3]*p[4]

println(newx)
println(newy)
println(newz)
