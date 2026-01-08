function trtrsRow(L, x) # This is the slow implementation
    n = length(x)
    for i=1:n
        for j = 1:i-1
            x[i] -= L[i,j] * x[j]
        end
        x[i] = x[i] / L[i,i]
    end
end

function trtrs(L, x)
    n = length(x)
    for j = 1:n
        x[j] = x[j] / L[j,j]
        for i = j+1:n
            x[i] -= L[i,j] * x[j]
            # This is fast because the row index is the inner loop
        end
    end
end
