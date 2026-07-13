# Working around my inability to come up with algorithms
# A collection of functions returning sequences of flops for 
# generating all binary functios at the given Hamming distance

# USAGE:
#
#   hamseq[depth](length)
#
# where depth is the Hamming distance and length is the number
# of bits in the string
# OUTPUT:
#   [depth, C] array of flip indices, where C is the number of
#   strings at the given HD

# NOTE: Currently `depth` is limited by 5
# TODO: make a universal version (only few usages, none for serious)


function hamseq1(NB)
    out = [i for i in 1:NB]
    return Int.(out)
end

function hamseq2(NB)
	NF = 2
	conf = zeros(NF,1)
	out = zeros(NF, 1)

	for conf[1] in 1:(NB - 1)
		Nk = 0
		for conf[2] in (conf[1]+1):(NB - Nk)
			out = [out conf]
		end
	end
    return Int.(out[:,2:end])
end

function hamseq3(NB)
	NF = 3
	conf = zeros(NF,1)
	out = zeros(NF, 1)

	for conf[1] in 1:(NB - NF + 1)
		Nk = NF - 2
	
		for conf[2] in (conf[1]+1):(NB - Nk)
			Nk = NF - 3

			for conf[3] in (conf[2]+1):(NB - Nk)
				out = [out conf]
			end
    
		end
	end
    return Int.(out[:,2:end])
end

function hamseq4(NB)
	NF = 4
	conf = zeros(NF,1)
	out = zeros(NF, 1)

	for conf[1] in 1:(NB - NF + 1)
		Nk = NF - 2
	
		for conf[2] in (conf[1]+1):(NB - Nk)
			Nk = NF - 3

			for conf[3] in (conf[2]+1):(NB - Nk)
				Nk = NF - 4

				for conf[4] in (conf[3]+1):(NB - Nk)
					out = [out conf]
				end
			end    
		end
	end
    return Int.(out[:,2:end])
end

function hamseq5(NB)
	NF = 5
	conf = zeros(NF,1)
	out = zeros(NF, 1)

	for conf[1] in 1:(NB - NF + 1)
		Nk = NF - 2
	
		for conf[2] in (conf[1]+1):(NB - Nk)
			Nk = NF - 3

			for conf[3] in (conf[2]+1):(NB - Nk)
				Nk = NF - 4

				for conf[4] in (conf[3]+1):(NB - Nk)
					Nk = NF - 5

					for conf[5] in (conf[4]+1):(NB - Nk)
						out = [out conf]
    				end
    			end
    		end
		end
	end
    return Int.(out[:,2:end])
end

hamseq = [hamseq1, hamseq2, hamseq3, hamseq4, hamseq5]
