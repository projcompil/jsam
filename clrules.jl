import OpenCL
const cl = OpenCL

const sum_kernel = "
   __kernel void sum(__global const float *a,
                     __global const float *b, 
                     __global float *c)
    {
      int gid = get_global_id(0);
      c[gid] = a[gid] + b[gid];
    }
"
a = rand(Float32, 50_000)
b = rand(Float32, 50_000)

device, ctx, queue = cl.create_compute_context()

a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=b)
c_buff = cl.Buffer(Float32, ctx, :w, length(a))

p = cl.Program(ctx, source=sum_kernel) |> cl.build!
k = cl.Kernel(p, "sum")

cl.call(queue, k, size(a), nothing, a_buff, b_buff, c_buff)

r = cl.read(queue, c_buff)

if isapprox(norm(r - (a+b)), zero(Float32))
    info("Success!")
else
    error("Norm should be 0.0f")
end

const max_kernel = "
   __kernel void sum(__global const int *prod,
   					 __global const int l,
                     __global float *input)
    {
		int gid = get_global_id(0);
		int maxi j;
		if(gid % l == 0) {
		maxi= prod[gid * l];
		for(j = gid ; j < gid + l ; j++) {
			maxi = max(maxi, prod[j]);
		}
		for(j = gid ; j < gid + l ; j++) {
			if(prod[j] == maxi) {
				input[j] = ;
			}
			else {
				input[j] = ;
			}
		}
		}
    }
"
const prod_kernel = "
   __kernel void sum(__global const int gamma,
                     __global const uint8 *b, 
                     __global int *c)
    {
      int gid = get_global_id(0);
      c[gid] = gamma * b[gid];
    }
"
const prod_kernel = "
   __kernel void sum(__global const int gamma,
                     __global const uint8 *b, 
                     __global int *c)
    {
      int gid = get_global_id(0);
      c[gid] = gamma * b[gid];
    }
"
