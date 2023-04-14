
list = [1, 2, 3, 4, 5, 6]
println(list)

for i in range(1, length(list))
    if list[i] == 2
        list[i] = 1
    end
end

println(list)