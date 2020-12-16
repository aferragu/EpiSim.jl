#utilities

function compute_harvard_levels(population::Number)

    pop = population/100000
    levels = (1,10,25,45)

    return levels.*pop

end
