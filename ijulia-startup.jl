using Pkg
Pkg.add("WebIO")
Pkg.build("WebIO")
using WebIO, IJulia
WebIO.install_jupyter_labextension()

# Run 'include("ijulia-startup.jl")' in julia, to install all the things that are needed.
