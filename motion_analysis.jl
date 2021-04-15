using Printf

# ------------------------------------
# output results
# ------------------------------------
function output_result(t, x, np, outdir)
    
    t = string(t)
    while length(t) < 8
        t = "0"*t
    end
    
    fff = outdir * "/result_" *t* ".dat"
    open(fff,"w") do f
        write(f,"result:x, y, z\n")
        for i in 1:np
            for l in 1:3
                a = @sprintf("%8.8e", x[i,l])
                write(f, a*" ")
            end
            write(f, "\n")
        end
    end
    println("\nwrite "*fff)
end

# ------------------------------
#  make directory
# ------------------------------
function make_dir(outdir)
    k = 0
    try rm(outdir,recursive=true)
    catch
        mkdir(outdir)
        k = 1
    end

    if k == 0
        mkdir(outdir)
    end
end

# ------------------------------
#  flux
# ------------------------------
function set_flux(Qbase, flux, m, Ix, Iy, Iz, Fx, Fy, Fz, Mx, My, Mz)

    flux[1] = Fx / m
    flux[2] = Fy / m
    flux[3] = Fz / m
    flux[4] = Qbase[1]  # v_x
    flux[5] = Qbase[2]  # v_y
    flux[6] = Qbase[3]  # v_z
    flux[7] = ( Iy - Iz ) * Qbase[8] * Qbase[9] / Ix + Mx / Ix # (Iy - Iz) omega_y omega_z / Ix + Mx/Ix
    flux[8] = ( Iz - Ix ) * Qbase[9] * Qbase[7] / Iy + My / Iy # (Iz - Ix) omega_z omega_x / Iy + My/Iy
    flux[9] = ( Ix - Iy ) * Qbase[7] * Qbase[8] / Iz + Mz / Iz # (Ix - Iy) omega_x omega_y / Iz + Mz/Iz
    #flux[10] = 0.5*(Qbase[13]*Qbase[7] - Qbase[12]*Qbase[8] + Qbase[11]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    #flux[11] = 0.5*(Qbase[12]*Qbase[7] + Qbase[13]*Qbase[8] - Qbase[10]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    #flux[12] = 0.5*(-Qbase[11]*Qbase[7] + Qbase[10]*Qbase[8] + Qbase[13]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    #flux[13] = 0.5*(-Qbase[10]*Qbase[7] - Qbase[11]*Qbase[8] - Qbase[12]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2

    #q0 = 0.5*( Qbase[13]*Qbase[7] - Qbase[12]*Qbase[8] + Qbase[11]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    #q1 = 0.5*( Qbase[12]*Qbase[7] + Qbase[13]*Qbase[8] - Qbase[10]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    #q2 = 0.5*(-Qbase[11]*Qbase[7] + Qbase[10]*Qbase[8] + Qbase[13]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    #q3 = 0.5*(-Qbase[10]*Qbase[7] - Qbase[11]*Qbase[8] - Qbase[12]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2

    q0 = 0.5*( Qbase[13]*Qbase[7] + Qbase[12]*Qbase[8] - Qbase[11]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    q1 = 0.5*(-Qbase[12]*Qbase[7] + Qbase[13]*Qbase[8] + Qbase[10]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    q2 = 0.5*( Qbase[11]*Qbase[7] - Qbase[10]*Qbase[8] + Qbase[13]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2
    q3 = 0.5*(-Qbase[10]*Qbase[7] - Qbase[11]*Qbase[8] - Qbase[12]*Qbase[9])  # (kai*omega_x - zeta*omega_y + eta*omega_z) / 2

    n = ( q0^2 + q1^2 + q2^2 + q3^2 )^0.5 

    flux[10] = q0
    flux[11] = q1
    flux[12] = q2
    flux[13] = q3
    #=
    if n == 0.0
        flux[10] = 0.0
        flux[11] = 0.0
        flux[12] = 0.0
        flux[13] = 0.0
    else
        flux[10] = q0 / n
        flux[11] = q1 / n
        flux[12] = q2 / n
        flux[13] = q3 / n
    end
    =#
    return flux
end

# ------------------------------
#  rotation
# ------------------------------
function rotation(x,y,z, Qbase)
    A = zeros(3, 3)
    q0 = Qbase[10]
    q1 = Qbase[11]
    q2 = Qbase[12]
    q3 = Qbase[13]

    #=

    A[1,1] = q0^2 + q1^2 - q2^2 - q3^2
    A[1,2] = 2*(q1*q2 + q0*q3)
    A[1,3] = 2*(q1*q3 - q0*q2)
    A[2,1] = 2*(q1*q2 - q0*q3)
    A[2,2] = q0^2 - q1^2 + q2^2 - q3^2
    A[2,3] = 2*(q2*q3 + q0*q1)
    A[3,1] = 2*(q1*q3 + q0*q2)
    A[3,2] = 2*(q2*q3 - q0*q1)
    A[3,3] = q0^2 - q1^2 - q2^2 + q3^2
    =#

    q1 = Qbase[13]
    q2 = Qbase[10]
    q3 = Qbase[11]
    q4 = Qbase[12]


    A[1,1] = 1 - 2*q3^2 - 2*q4^2
    A[1,2] = 2*q2*q3 + 2*q1*q4
    A[1,3] = 2*q2*q4 - 2*q1*q3

    A[2,1] = 2*q2*q3 - 2*q1*q4
    A[2,2] = 1 - 2*q4^2 - 2*q2^2
    A[2,3] = 2*q3*q4 + 2*q1*q2

    A[3,1] = 2*q2*q4 + 2*q1*q3
    A[3,2] = 2*q3*q4 - 2*q1*q2
    A[3,3] = 1 - 2*q2^2 - 2*q3^2


    newx = A[1,1]*x + A[1,2]*y + A[1,3]*z
    newy = A[2,1]*x + A[2,2]*y + A[2,3]*z
    newz = A[3,1]*x + A[3,2]*y + A[3,3]*z    
    return newx, newy, newz
end

# ------------------------------
#  main 
# ------------------------------
function main()

    m = 1.0 # [kg]

    Fx = 0.0 # [N]
    Fy = 0.0 # [N]
    Fz = 0.0 # [N]
    Mx = 0.0 # [N/m]
    My = 0.0 # [N/m]
    Mz = 0.0e0 # [N/m]

    init_v_x = 0.0
    init_v_y = 0.0
    init_v_z = 0.0
    init_gc_x = 0.0
    init_gc_y = 0.0
    init_gc_z = 0.0
    init_omega_x = 1.5e-4
    init_omega_y = 0.0
    init_omega_z = 0.0

    dt = 1.0e-6 # [s]
    nt = 1.0e6
    outstep = 1.0e4

    outdir = "post_result"
    make_dir(outdir)
    outdir = "result"
    make_dir(outdir)

    # quontanion 
    # q = (1, 0, 0, 0)
    # q = (ux sin(theta/2), uy sin(theta/2), uz sin(theta/2), cos(theta/2))
    # これにより初期角度と軸を決定できる
    # u = (ux, uy, uz) ：回転軸方向の単位ベクトル
    init_q_0 = 0.0
    init_q_1 = 0.0
    init_q_2 = 0.0
    init_q_3 = 1.0

    
    # ------------------------------
    #  make block
    # ------------------------------
    np = 8
    x = zeros(np, 3)

    for i in 1:np
        if (i-1) % 8 < 4
            x[i, 1] = 1.0
        else
            x[i, 1] = -1.0
        end

        if (i-1) % 4 < 2
            x[i, 2] = 1.0
        else
            x[i, 2] = -1.0
        end

        if (i-1) % 2 < 1
            x[i, 3] = 1.0
        else
            x[i, 3] = -1.0
        end
    end
    
    # ------------------------------
    #  primitive
    # ------------------------------
    # vecVelo, GC, vecOmega, くっぉーたにオン
    nprim = 13
    Qbase = zeros(nprim)
    Qbase[1] = init_v_x
    Qbase[2] = init_v_y
    Qbase[3] = init_v_z
    Qbase[4] = init_gc_x
    Qbase[5] = init_gc_y
    Qbase[6] = init_gc_z
    Qbase[7] = init_omega_x
    Qbase[8] = init_omega_y
    Qbase[9] = init_omega_z
    Qbase[10] = init_q_0
    Qbase[11] = init_q_1
    Qbase[12] = init_q_2
    Qbase[13] = init_q_3

    flux = zeros(nprim)
    # ------------------------------
    # innertia moment
    # ------------------------------
    # 長さ2より
    Ix = m * (2^2 + 2^2) / 12
    Iy = m * (2^2 + 2^2) / 12
    Iz = m * (2^2 + 2^2) / 12

    # ------------------------------
    # main loop
    # ------------------------------
    for t in 1:Int(nt)

        flux = set_flux(Qbase, flux, m, Ix, Iy, Iz, Fx, Fy, Fz, Mx, My, Mz)
        
        # Euler explicit
        for i in 1:nprim
            Qbase[i] = Qbase[i] + dt*flux[i]           
        end

        # 規格化
        #=
        q0 = Qbase[10]
        q1 = Qbase[11]
        q2 = Qbase[12]
        q3 = Qbase[13]
        n = ( q0^2 + q1^2 + q2^2 + q3^2 )^0.5 
        Qbase[10] = q0*n
        Qbase[11] = q1*n
        Qbase[12] = q2*n
        Qbase[13] = q3*n
        =#

        #println(Qbase)


        # move all point
        for i in 1:np
            x[i,1] += Qbase[1] * dt
            x[i,2] += Qbase[2] * dt
            x[i,3] += Qbase[3] * dt

            # rotation
            x[i,1], x[i,2], x[i,3] = rotation(x[i,1], x[i,2], x[i,3], Qbase)
        end
        
        #=
        println("throw(UndefVarError(:x))")
        println(Qbase)
        println(flux)
        println(x)    
        =#

        if t % Int(outstep) == 0 
            output_result(t, x, np, outdir)
        end

        #=

        if t == 100
            Qbase[7] = -1.5e-3
        elseif t == 200
            Qbase[7] = 1.5e-3
        else
            Qbase[7] = 0.0
        end
        =#
        
        #=

        q0 = Qbase[10]
        q1 = Qbase[11]
        q2 = Qbase[12]
        q3 = Qbase[13]

        theta = acos(q3)*2
        ux    = q0 / sin(theta/2) 
        uy    = q1 / sin(theta/2) 
        uz    = q2 / sin(theta/2) 
        println(" ")
        println(theta)
        println(ux)
        println(uy)
        println(uz)
        #throw(UndefVarError(:x))
        =#


        
        if isequal(x[1,1], NaN) == true
            println("\n -------------------- \n ")
            println(" ")
            println(" diverge ")
            println(" ")
            println("\n -------------------- \n ")
            throw(UndefVarError(:x))
        end

        if abs(x[1,1]) > 2
            println("\n -------------------- \n ")
            println(" ")
            println(" use small dt ")
            println(" ")
            println("\n -------------------- \n ")
            throw(UndefVarError(:x))
        end

    end
    
end
# ------------------------------
# ------------------------------
main()