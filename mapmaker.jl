require("Winston")
require("Winston/src/Plot")
require("suitesparse")

module mapmaker
import Main.Plot, Main.Winston
export A,cost,control,mapinfo,main

type MapInfo
    fname::String
    west::Int
    east::Int
    north::Int
    south::Int
    W::Int
    H::Int
end
type MapPoint
    north::Float64
    east::Float64
end
function ntoh!(a::Array)
    for i = 1:numel(a)
        a[i] = ntoh(a[i])
    end
    a
end

function show_cost(cost)
    (W,H) = size(cost)
    p = Plot.imagesc((0,W),(H,0),cost)
    Plot.tk(p)
end
function show_control(control)
    (W,H) = size(control)
    p = Plot.imagesc((0,W),(H,0),control,(0,8))
    Plot.tk(p)
end
function show_map(A,mapinfo,dest)
    W = mapinfo.W
    H = mapinfo.H
    x0 = [int((dest.east - mapinfo.west) / (mapinfo.east - mapinfo.west) * W),
          int((dest.north - mapinfo.south) / (mapinfo.north - mapinfo.south) * H)] 
    p=Plot.imagesc((0,W),(H,0),A);
    add(p,Winston.Point(x0[1],x0[2],"symboltype","asterisk"));
    add(p,Winston.Point(x0[1],x0[2],"symboltype","circle"));
    Plot.tk(p,"width",1440,"height",720)
end


function policy_iter(x0,W,H,alpha,control,nextx,u_cost,g)
    T = speye(W*H,W*H)
    g[:] = 0
    for i = 1:W
        for j = 1:H
            ci = control[j,i]
            if ci != 0
                T[j+(i-1)*H, nextx[ci,j,i]] = -alpha
                g[j+(i-1)*H] = u_cost[ci,j,i]
            end
        end
    end
    J = T\g
    cost = reshape(J,H,W)
    cost
end

function value_iter(x0,W,H,alpha,prevcost,nextx,u_cost,control,cost)
    changed = false
    for i = 1:W
        for j = 1:H
            if x0[2] == j && x0[1] == i
                control[j,i] = 0
                cost[j,i] = 0
            else
                costs = u_cost[:,j,i] + alpha*prevcost[nextx[:,j,i]]
                (costi, ci) = findmin(costs)
                if !changed && control[j,i] != ci
                    changed = true
                end
                control[j,i] = ci
                cost[j,i] = costi
            end
        end
    end
    changed
end

## Policy Controls Key:
#    812
#    703
#    654
function main(A,mapinfo,dest,jump,alpha,piter_range,viter_range,vtol)
    A=A[1:jump:end,1:jump:end]
    W = int(mapinfo.W/jump)
    H = int(mapinfo.H/jump)
    x0 = [int((dest.east - mapinfo.west) / (mapinfo.east - mapinfo.west) * W),
          int((dest.north - mapinfo.south) / (mapinfo.north - mapinfo.south) * H)]
    println("main conditions: W: $W, H: $H")
    println("    jump: $jump, alpha: $alpha, vtol: $vtol")
    println("    piter: $piter_range, viter: $viter_range")

    
    ## Init grid (direct) and costs (euclidean)
    cost = zeros(H,W)
    control = zeros(Int,H,W)
    for i = 1:W # E/W
        for j = 1:H # N/S
            #ci = mod(int(atan2(x0[2]-j, x0[1]-i) / 2 / pi * 8 + 6),8)+1
            #costi = hypot(j-x0[2], i-x0[1])
            if     j <= x0[2] && i <  x0[1]
                ci = 7
            elseif j <  x0[2] && i >= x0[1]
                ci = 5
            elseif j >= x0[2] && i > x0[1]
                ci = 3
            elseif j >  x0[2] && i <= x0[1]
                ci = 1
            else
                ci = 0
            end
            costi = abs(j-x0[2]) + abs(i-x0[1])
            control[j,i] = ci
            cost[j,i] = costi
        end
    end
    i = x0[1]
    j = x0[2]
    control[j,i] = 0
    cost[j,i] = 0
    
    ## Show stats
    show_cost(cost)
    show_control(control)
    #return(control,cost)

    # Prep iteration
    multiplier = [1,sqrt(2),1,sqrt(2),1,sqrt(2),1,sqrt(2)]
    nextx = zeros(Int,8,H,W)
    u_cost = zeros(8,H,W)
    for i = 1:W
        for j = 1:H
            h = A[j,i]
            nx_xs = [
                mod(j-2,H)+1  i             # down
                mod(j-2,H)+1  mod(i-2,W)+1  # down left
                j             mod(i-2,W)+1  # left
                mod(j,H)+1    mod(i-2,W)+1  # up left
                mod(j,H)+1    i             # up
                mod(j,H)+1    mod(i,W)+1    # up right
                j             mod(i,W)+1    # right
                mod(j-2,H)+1  mod(i,W)+1    # down right
                ]
            nx_x = nx_xs[:,1]+(nx_xs[:,2]-1)*H
            nextx[:,j,i] = nx_x
            next_h = A[nx_x]
            # shortest path cost
            # u_cost[:,j,i] = multiplier

            # small linear cost fcn gives nice solution, approx to cost map
            # u_cost[:,j,i] = (abs((next_h-h)/800) + 1).*multiplier
            
            # medium square cost fcn gives nicer solution, approx flow field
             u_cost[:,j,i] = (((next_h-h)/1000).^2 + 1).*multiplier
        end
    end
   
    ## Policy/Value iteration of the actual, best policy/costs
    local piter, viter_sum = 0
    try
    g = zeros(H*W)
    for piter = piter_range
        local changed,viter
        if piter != 0
            @time begin
                cost = policy_iter(x0,W,H,alpha,control,nextx,u_cost,g)
                prev_cost = copy(cost)
                #show_cost(cost)
                show_control(control)
                print("piteration $piter ")
            end
        end
        @time begin
            for viter = viter_range
                @time begin
                    prev_cost = copy(cost)
                    changed = value_iter(x0,W,H,alpha,prev_cost,nextx,u_cost,control,cost)
                    #show_cost(cost)
                    show_control(control)
                    print("  viteration $viter ")
                end
                if !changed && norm(prev_cost - cost) < W*H*vtol
                    println("probably reached convergance $(norm(cost))")
                    break
                    #error(InterruptException())
                end
            end
            viter_sum += viter
            print("viter total ")
        end 
        if !changed && viter <= 1
            println("reached convergance $(norm(cost))")
            break
        end
    end
    catch e
        if !isa(e,InterruptException)
            rethrow(e)
        end
    end
    ## Show stats
    show_cost(cost)
    show_control(control)
    println("total value iterations: $viter_sum")
    println("total policy iterations: $piter")

    cost,control
end

## Load data
function read_map(mapinfo)
    f=open(mapinfo.fname)
    A = Array(Int16,mapinfo.W,mapinfo.H)
    read(f, A)
    ntoh!(A)
    close(f)
    A
end

# source: http://pds-geosciences.wustl.edu/missions/mgs/megdr.html
const mapinfo = MapInfo("megt90n000cb.img", 0,360, 90,-90, 1440,720)
const A = read_map(mapinfo)'

# source: http://toolserver.org/~geohack/geohack.php?pagename=Mars_Science_Laboratory&params=4.5895_S_137.4417_E_globe:Mars
const curiosity = MapPoint(-4.5895,137.4417)

function main(jump,alpha,piter,viter,vtol)
    main(A, mapinfo, curiosity, jump, alpha, piter, viter, vtol)
end

function main()
    (cost,control) = @time main(16, 1, 1:1000, 1:1, 1e9) # (optimistic) policy iteration
    #(cost,control) = @time main(16, 1, 0, 1:10000, 1e-9) # value iteration
end

function mapkey()
    (cost,control) = main(zeros(9,9), MapInfo("", 1,3, 1,3, 9,9), MapPoint(2.,2.), 1, 1, 1:0, 1:0, 1e-9)
end

end

