import ArgParse

function input(prompt="")
	   print(prompt)
	   chomp(readline())
end

function parse_commandline()
	   table = ArgParse.ArgParseSettings()
	   ArgParse.@add_arg_table table begin
			   "--fsum", "-f"
			   help = "Function that will be used to retrieve corrupted inputs. 0 for sum_of_sum!. 1 for sum_of_max!"
			   arg_type = Int
			   default = 0

			   "-m"
			   help = "Number of messages stored in the network."
			   arg_type = Int
			   default = 5000
			   
			   "-l"
			   help = "Number of neurons per cluster."
			   arg_type = Int 
			   default = 128

			   "-c"
			   help = "Number of clusters."
			   arg_type = Int 
			   default = 8

			   "--gamma", "-g"
			   help = "Memory parameter."
			   arg_type = Int 
			   default = 1

			   "--erasures", "-e"
			   help = "Number of erasures."
			   arg_type = Int 
			   default = 4

			   "--proba-cons", "-p"
			   help = "Probability of adding connection."
			   arg_type = Float64
			   default = 1.0

			   "--proba-des"
			   help = "Probability of dropping connection."
			   arg_type = Float64
			   default = 0.0

			   "--degree", "-d"
			   help = "Number of edges per vertex. < 1 means all."
			   arg_type = Int 
			   default = 0

			   "--activities", "-a"
			   help = "Number of activated neurons per cluster."
			   arg_type = Int 
			   default = 1

			   "--diffusion"
			   help = "Number of diffusion steps before contraction."
			   arg_type = Int 
			   default = 1

			   "--iterations", "-i"
			   help = "Maximum number of iterations to retrieve a message."
			   arg_type = Int 
			   default = 4

			   "--tests", "-t"
			   help = "Number of tests to assess the score."
			   arg_type = Int 
			   default = 1000

			   "--trials", "-n"
			   help = "Number of trials simulated"
			   arg_type = Int 
			   default = 1

			   "--fcorrupt", "-z"
			   help = "Cluster corruption function. 0 for erasing clusters. 1 for shuffling them."
			   arg_type = Int 
			   default = 0

			   "--winners", "-w"
			   help = "Number of winners per cluster."
			   arg_type = Int 
			   default = 1

			   "--nowrite"
			   help = "Do not write results in file."
			   action = :store_true

			   "--scan"
			   help = "Test all erasures or errors from 0 to c."
			   action = :store_true

			   "--lscan"
			   help = "Scan on sizes of clusters. (hardcoded at the moment)."
			   action = :store_true

			   "--mscan"
			   help = "Scan on numbers of messages (hardcoded at the moment)."
			   action = :store_true

			   "--file"
			   help = "Prefix of the output file."
			   default = "clique_results"

			   "--fconfig"
			   help = "Name of configuration file."
			   default = "config.json"

			   "--config", "-C"
			   help = "Use a configuration file."
			   action = :store_true

			   "--dir"
			   help = "Directory of the output file.."
			   default = "results"

			   "--quiet", "-q"
			   help = "Outputs less text. (Not implemented at the moment)."
			   action = :store_true

			   "--bind-to"
			   help = "Accept this option."
			   arg_type = String
			   "--worker"
			   help = "Accept this option."
			   action = :store_true

	   end
	   return ArgParse.parse_args(table)
end
