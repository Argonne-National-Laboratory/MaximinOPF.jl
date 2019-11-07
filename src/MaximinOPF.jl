module MaximinOPF
using PowerModels

greet() = print("Hello World!")

function MaximinOPFModel(case, powerform, nLine)
	println("Hello MaximinOPFModel")
	println(case)
	println(powerform)
	println(nLine)

	pm = build_model(case, powerform, PowerModels.post_opf)

	return pm
end

end # module
