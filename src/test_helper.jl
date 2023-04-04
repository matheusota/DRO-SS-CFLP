include("data.jl")
include("file_handler.jl")
include("solution.jl")
include("distributions_helper.jl")

const DOCKER_CMD_1 = `docker run --rm -v /home/matheus/Projects/StochasticVRP/SVRP/CCVRP_Julia:/CVRP bapdock /CVRP/src/run.jl `
const DOCKER_CMD_2 = " -m 5 -M 5 -u 950 -o "
const TMP_INPUT = "tmp/input/tmp_input.vrp"
const TMP_OUTPUT = "tmp/output/tmp_output"

const DOCKER_ARGS = ["run",
"--rm",
"-v",
"/home/matheus/Projects/StochasticVRP/SVRP/CCVRP_Julia:/CVRP",
"bapdock",
"/CVRP/src/run.jl",
"/CVRP/" * TMP_INPUT,
"-m",
"5",
"-M",
"5",
"-u",
"950",
"-o",
"/CVRP/" * TMP_OUTPUT * ".sol",
"-t",
"/CVRP/" * TMP_OUTPUT * ".tex"]

const OUT_OF_SAMPLES = 10000

function failures_from_SVRP_Docker(distributionName::String, varianceType::String, fileString::String)
	data = read_instance("", fileString)

	# sample scenarios from the chosen distribution
	if distributionName == "Independent Normal"
		distribution = get_normal_independent_demands(data, varianceType, "")
	elseif distributionName == "Joint Normal"
		distribution = get_normal_joint_demands(data, varianceType, "")
	end

	data.scenarios = get_scenarios_from_distribution(distribution, data.n, 20)
	data.nScenarios = 20

	# write to temporary input file
	write_instance(data, TMP_INPUT)

	# run algorithm and create output image
	run(`docker $(DOCKER_ARGS)`)
	run(`pdflatex $([TMP_OUTPUT * ".tex", "-output-directory", "tmp/output"])`)
	run(`pdftoppm $(["-jpeg", "-r", "100", TMP_OUTPUT * ".pdf", TMP_OUTPUT])`)
	println("FINISHED RUNNING")

	# get number of failures
	solution = readsolution(Dict{String,Any}("sol" => TMP_OUTPUT * ".sol", "round" => false))
	inSampleFailures = test_for_failures(data, solution, distribution, true)
	outOfSampleFailures = test_for_failures(data, solution, distribution, false)

	write_failures(inSampleFailures, outOfSampleFailures, TMP_OUTPUT * "_failures.txt")
end

function run_test()
	failures_from_SVRP_solution("./data/A/A-n32-k5-scen-L666.svrp", "./out/A-n32-k5-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n32-k5-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/A/A-n34-k5-scen-L666.svrp", "./out/A-n34-k5-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n34-k5-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/A/A-n36-k5-scen-L666.svrp", "./out/A-n36-k5-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n36-k5-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/A/A-n37-k5-scen-L666.svrp", "./out/A-n37-k5-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n37-k5-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/B/B-n39-k5-scen-L666.svrp", "./out/B-n39-k5-scen-L666.sol", "Joint Normal", "Low", "./data/B/B-n39-k5-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/A/A-n44-k6-scen-L666.svrp", "./out/A-n44-k6-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n44-k6-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/A/A-n45-k6-scen-L666.svrp", "./out/A-n45-k6-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n45-k6-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/P/P-n50-k10-scen-L666.svrp", "./out/P-n50-k10-scen-L666.sol", "Joint Normal", "Low", "./data/P/P-n50-k10-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/P/P-n51-k10-scen-L666.svrp", "./out/P-n51-k10-scen-L666.sol", "Joint Normal", "Low", "./data/P/P-n51-k10-jnorm-L.svrp")
	failures_from_SVRP_solution("./data/A/A-n55-k9-scen-L666.svrp", "./out/A-n55-k9-scen-L666.sol", "Joint Normal", "Low", "./data/A/A-n55-k9-jnorm-L.svrp")
end

function failures_from_SVRP_solution(inputPath::String, outputPath::String, distributionName::String, varianceType::String, distributionFile::String)
	if length(inputPath) > 50
		data = read_instance("", inputPath)
	else
		data = read_instance(inputPath)
	end

	# sample scenarios from the chosen distribution
	if distributionName == "Independent Normal"
		distribution = get_normal_independent_demands(data, varianceType, distributionFile)
	elseif distributionName == "Joint Normal"
		distribution = get_normal_joint_demands(data, varianceType, distributionFile)
	end

	solution = readsolution(Dict{String,Any}("sol" => outputPath, "round" => false))
	inSampleFailures, inSampleCost = test_for_failures(data, solution, distribution, true)
	outOfSampleFailures, outOfSampleCost = test_for_failures(data, solution, distribution, false)

	println("Input File: $(inputPath)")
	println("In Sample Failures: $(inSampleFailures) / ($(data.k) * $(data.nScenarios)) = $(inSampleFailures / (data.k * data.nScenarios))")
	println("In Sample Extra Cost: $(inSampleCost)")
	println("Out of Sample Failures: $(outOfSampleFailures) / ($(data.k) * $(OUT_OF_SAMPLES)) = $(outOfSampleFailures / (data.k * OUT_OF_SAMPLES))")
	println("Out of Sample Extra Cost: $(outOfSampleCost)")

	open("out/sample_table.txt", "a") do f
		write(f, "$inputPath & $(inSampleFailures / (data.k * data.nScenarios)) & $(inSampleCost) & $(outOfSampleFailures / (data.k * OUT_OF_SAMPLES)) & $(outOfSampleCost) \\\\ \n")
	end
end

function test_for_failures(data::DataCVRP, solution::Solution, distribution, inSample)
	# get scenarios to be tested
	if inSample
		scenarios = data.scenarios
		nScenarios = data.nScenarios
	else
		scenarios = get_scenarios_from_distribution(distribution, data.n, OUT_OF_SAMPLES)
		nScenarios = OUT_OF_SAMPLES
	end

	# simulate executing the routes to get the failures
	countFailures = 0
	extraCost = 0.0
	for i in 1:nScenarios
		scenario = zeros(Float64, data.n)
		for j in 1:data.n
			scenario[j] = scenarios[j][i]
		end

		nFailures, cost = get_failures(data, solution, data.Q, scenario)
		countFailures += nFailures
		extraCost += cost
	end

	return countFailures, extraCost / nScenarios
end

function write_failures(inSampleFailures::Int64, outOfSampleFailures::Int64, filePath::String)
	open(filePath, "w") do f
		write(f, "In: $(inSampleFailures)\n")
		write(f, "Out: $(outOfSampleFailures)\n")
	end
end

function get_failures(data::DataCVRP, solution::Solution, capacity::Float64, demands::Array{Float64})
	nFailures = 0
	extraCost = 0.0

	for route in solution.routes
		s = ""
		s *= "route: $route\n"
		curr = 0.0
		for i in route
			curr += demands[i + 1]
			s *= "$(demands[i + 1])\n"
			if curr > capacity
				s *= "curr: $curr \n"
				nFailures += 1
				extraCost += distance(data, (0, i))
				curr -= capacity
				#println(s)
			end
		end
	end

	return nFailures, extraCost
 end

 run_test()
