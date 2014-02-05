

%sum of surrounding voxels
binsum = zeros(size(binedges));
counter = 0;
for z=2:(zdim-1)
    for i=2:(xdim-1)
        for j=2:(ydim-1)
            if binedges(j,i,z) == 1
            for zz=-1:1 
                for ii=-1:1
                    for jj=-1:1
                        counter = counter + binedges((j+jj),(i+ii),(z+zz));
                    end
                end
            end
            binsum(j,i,z) = 27-counter; %9 to 27 if 3d
            counter = 0;

            end
        end
    end
end

%compare to surrounding binsum, mark local max (so less filled in nearby
%voxels, so edge of a hole?)

localmax = binedges;

for z=2:(zdim-1)
    for i=2:(xdim-1)
        for j=2:(ydim-1)
            if binedges(j,i,z) == 1 &&...
                    binsum(j,i,z)>binsum((j-1),(i-1),(z)) &&...
                    binsum(j,i,z)>binsum((j-1),(i),(z)) &&...
                    binsum(j,i,z)>binsum((j-1),(i+1),(z)) &&...
                    binsum(j,i,z)>binsum((j),(i-1),(z)) &&...
                    binsum(j,i,z)>binsum((j),(i+1),(z)) &&...
                    binsum(j,i,z)>binsum((j+1),(i-1),(z)) &&...
                    binsum(j,i,z)>binsum((j+1),(i),(z)) &&...
                    binsum(j,i,z)>binsum((j+1),(i+1),(z)) &&...
                    binsum(j,i,z)>binsum((j-1),(i-1),(z-1)) &&...
                    binsum(j,i,z)>binsum((j-1),(i),(z-1)) &&...
                    binsum(j,i,z)>binsum((j-1),(i+1),(z-1)) &&...
                    binsum(j,i,z)>binsum((j),(i-1),(z-1)) &&...
                    binsum(j,i,z)>binsum((j),(i),(z-1)) &&...
                    binsum(j,i,z)>binsum((j),(i+1),(z-1)) &&...
                    binsum(j,i,z)>binsum((j+1),(i-1),(z-1)) &&...
                    binsum(j,i,z)>binsum((j+1),(i),(z-1)) &&...
                    binsum(j,i,z)>binsum((j+1),(i+1),(z-1)) &&...
                    binsum(j,i,z)>binsum((j-1),(i-1),(z+1)) &&...
                    binsum(j,i,z)>binsum((j-1),(i),(z+1)) &&...
                    binsum(j,i,z)>binsum((j-1),(i+1),(z+1)) &&...
                    binsum(j,i,z)>binsum((j),(i-1),(z+1)) &&...
                    binsum(j,i,z)>binsum((j),(i),(z+1)) &&...
                    binsum(j,i,z)>binsum((j),(i+1),(z+1)) &&...
                    binsum(j,i,z)>binsum((j+1),(i-1),(z+1)) &&...
                    binsum(j,i,z)>binsum((j+1),(i),(z+1)) &&...
                    binsum(j,i,z)>binsum((j+1),(i+1),(z+1))
                localmax(j,i,z) = 2;

             
            end
        end
    end
end