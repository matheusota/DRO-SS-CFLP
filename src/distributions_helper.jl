import Unicode
using Random, Distributions, LinearAlgebra

Random.seed!(666)

function get_deviations(data::DataCVRP, varianceType::String)
   if varianceType == "Low"
      a, b = 0.07, 0.13
   else
      a, b = 0.14, 0.26
   end

   return [rand(Uniform(a * v.demand, b * v.demand)) for v in data.G′.V′[2:end]]
end

function get_normal_independent_demands(data::DataCVRP, type::String, distributionFile::String)
   if distributionFile == ""
      deviations = get_deviations(data, type)
   else
      dataAux = read_instance(distributionFile)
      deviations = [sqrt(dataAux.variance[i]) for i in 2:data.n]
   end

   return [Normal(v.demand, deviations[v.id_vertex]) for v in data.G′.V′[2:end]]
end

function get_normal_joint_demands(data::DataCVRP, varianceType::String, distributionFile::String)
   if distributionFile == ""
      deviations = get_deviations(data, varianceType)
      gamma = [i != j ? 1.0 / (distance(data, (i, j)) * rand(Uniform(0.4, 1.6))) : 1.0 for i in 1:data.n-1, j in 1:data.n-1]
      max_gamma, min_gamma = maximum(gamma), minimum(gamma)
      correlations = [i <= j ? (0.2 * gamma[i, j]) / (max_gamma + min_gamma) : (0.2 * gamma[j, i]) / (max_gamma + min_gamma) for i in 1:data.n-1, j in 1:data.n-1]
      covariances = [correlations[i, j] * deviations[i] * deviations[j] for i in 1:data.n-1, j in 1:data.n-1]
   else
      dataAux = read_instance(distributionFile)
      covariances = dataAux.covariance
   end

   #display(covariances)
   means = [v.demand for v in data.G′.V′[2:end]]
   return [MvNormal(means, covariances)]
end

function get_scenarios_from_distribution(distribution, nVertices, nScenarios::Int)
   scenarios = [zeros(Cint, nScenarios) for i in 1:nVertices]

   for j in 1:nScenarios
      if length(distribution) == 1
         scenario_demands = rand(distribution[1])
      else
         scenario_demands = [rand(d) for d in distribution]
      end

      scenario_demands = [0; [round(d) for d in scenario_demands]]

      for i in 1:nVertices
         scenarios[i][j] = Int32(scenario_demands[i])
      end
   end

   return scenarios
end
