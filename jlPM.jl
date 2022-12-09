#This code solves the PM model subject to a given forcing with/without vorticity
using Plots
using DifferentialEquations

function PM(r,Mz,F)
    return -r*Mz + F
end

function solve_PM(t,Mz0,r,F,solver;ζ=zeros(length(t)))
    dt = t[2]-t[1]
    Mz = zeros(length(t)).+ 1im*zeros(length(t))
    Mz[1] = Mz[0]

    if solver = "FE"
        for i = 2:length(Mz)
            r_eff = r+1im*ζ[i-1]/2
            Mz[i] = Mz[i-1]+dt*PM(r_eff,Mz[i-1],F[i-1])
        end
        return Mz
    elseif solver = "CN"
        for i = 2:length(Mz)
            r_eff1 = r+1im*ζ[i-1]/2
            r_eff2 = r+1im*ζ[i]/2
            Mz[i] = ((1-r_eff1*dt/2)*Mz[i-1]+dt*(F[i]+F[i-1])/2)/(1+r_eff2*dt/2)
        end
        return Mz
    else 
        return "Not a valid solver type"
    end
end