function writecsvstep(t, x, ρ, u, v, p)
    file = open(joinpath("data", "data_$t.csv"), "w")
    write(file, "position,density,energy,velocity,pressure\n")
    writecsv(file, [x ρ u v p])
    close(file)
end
