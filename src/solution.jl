import Unicode

mutable struct Solution
   cost::Union{Int,Float64}
   routes::Array{Array{Int}}
end

# build Solution from the variables x
function getsolution(data::DataCVRP, x_values, objval, app::Dict{String,Any})
   E, dim = edges(data), dimension(data)
   adj_list = [[] for i in 1:dim] 
   for e in E
      val = x_values[e]
      if val > 0.5
         push!(adj_list[e[1]+1], e[2])
         push!(adj_list[e[2]+1], e[1])
         if val > 1.5
            push!(adj_list[e[1]+1], e[2])
            push!(adj_list[e[2]+1], e[1])
         end
      end
   end
   visited, routes = Dict([(i,false) for i in customers(data)]), []
   for i in adj_list[1] # starting from adjacent list of depot 0
      if !visited[i]
         r, prev = [0], 0
         push!(r, i)
         visited[i] = true
         length(adj_list[i+1]) != 2 && error("Problem trying to recover the route from the x values. "*
                                             "Customer $i has $(length(adj_list[i+1])) incident edges.")
         next, prev = (adj_list[i+1][1] == prev) ? adj_list[i+1][2] : adj_list[i+1][1], i
         push!(r, next)
         maxit, it = dim, 0
         while next != 0 && it < maxit # while depot 0 is not reached
            length(adj_list[next+1]) != 2 && error("Problem trying to recover the route from the x values. "* 
                                                   "Customer $next has $(length(adj_list[next+1])) incident edges.")
            visited[next] = true
            aux = next
            next, prev = (adj_list[next+1][1] == prev) ? adj_list[next+1][2] : adj_list[next+1][1], aux
            push!(r, next)
            it += 1
         end
         (it == maxit) && error("Problem trying to recover the route from the x values. "*
                                "Some route can not be recovered because the return to depot is never reached")
         push!(routes, r)
      end
   end 
   !isempty(filter(a->a==false,collect(values(visited)))) && error("Problem trying to recover the route from the x values. "*
                              "At least one customer was not visited or there are subtours in the solution x.")
   if app["round"]
      objval = trunc(Int, round(objval))
   end
   return Solution(objval, routes)
end

function print_routes(solution)
   for (i,r) in enumerate(solution.routes)
      print("Route #$i: ") 
      for j in r
         print("$j ")
      end
      println()
   end
end

# checks the feasiblity of a solution
function checksolution(data::DataCVRP, solution)
   dim, Q = dimension(data), veh_capacity(data)
   visits = Dict([(i,0) for i in customers(data)])
   sum_cost = 0.0
   for (i,r) in enumerate(solution.routes)
      sum_demand, prev = 0.0, r[1]
      for j in r[2:end-1]
         visits[j] += 1
         (visits[j] == 2) && error("Customer $j was visited more than once")
         sum_cost += distance(data, (prev,j))
         sum_demand += d(data, j)
         prev = j
      end
      sum_cost += distance(data, (prev,r[end]))
      (sum_demand > Q) && error("Route #$i is violating the capacity constraint. Sum of the demands is $(sum_demand) and Q is $Q")
   end
   !isempty(filter(a->a==0,collect(values(visits)))) && error("The following customers were not visited: $([k for (k,v) in visits if v==0])")
   (abs(solution.cost-sum_cost) > 0.001) && error("Cost calculated from the routes ($sum_cost) is different from that passed as"*
                                                                                                  " argument ($(solution.cost)).") 
end

# read solution from file (CVRPLIB format)
function readsolution(app::Dict{String,Any})
   if length(app["sol"]) > 50
      str = app["sol"]
   else
      str = read(app["sol"], String)
   end

   breaks_in = [' '; ':'; '\n';'\t';'\r']
   aux = split(str, breaks_in; limit=0, keepempty=false) 
   sol = Solution(0, [])
   j = 3
   while j <= length(aux)
      r = []
      while j <= length(aux)
         push!(r, parse(Int, aux[j]))
         j += 1
         if contains(lowercase(aux[j]), "cost") || contains(lowercase(aux[j]), "route")
            break
         end
      end
      push!(sol.routes, r)
      if contains(lowercase(aux[j]), "cost")
         if !app["round"]
            sol.cost = parse(Float64, aux[j+1])
         else 
            sol.cost = parse(Int, aux[j+1])
         end
         return sol
      end
      j += 2 # skip "Route" and "#j:" elements
   end
   error("The solution file was not read successfully. The file must be in the CVRPLIB format.")
   return sol
end

# write solution in a file
function writesolution(solpath, solution)
   open(solpath, "w") do f
      for (i,r) in enumerate(solution.routes)
         write(f, "Route #$i: ")
         for j in r
            write(f, "$j ") 
         end
         write(f, "\n")
      end
      write(f, "Cost $(solution.cost)\n")
   end
end

# write solution as TikZ figure (.tex) 
function drawsolution(tikzpath, data, solution)
   open(tikzpath, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
      # get limits to draw
      pos_x_vals = [i.pos_x for i in data.G′.V′]
      pos_y_vals = [i.pos_y for i in data.G′.V′]
      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")
      for cid in customers(data)
         i = data.G′.V′[cid+1]
         x_plot = scale_fac*i.pos_x
         y_plot = scale_fac*i.pos_y
         write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
      end
      # draw depot   
      i = data.G′.V′[1]
      x_plot = scale_fac*i.pos_x
      y_plot = scale_fac*i.pos_y
      write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=yellow, inner sep=0.05cm, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
      for r in solution.routes
         prev = r[1]
         for i in r[2:end-1]
            e = (prev,i)
            write(f, "\t\\draw[-,line width=0.8pt] (v$(e[1])) -- (v$(e[2]));\n")
            prev = i
         end
         write(f, "\t\\draw[-,line width=0.8pt] (v$(r[end])) -- (v$(prev));\n") 
      end
      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end   
end