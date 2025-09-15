function filename = append_flags(filename, crit1, crit2, crit3)
    if not(crit1 == 0)
        filename = filename + " - badMarker";
    end
    if not(crit2 == 0)
        filename = filename + " - badForcePlate";
    end
    if not(crit3 == 0)
        filename = filename + " - badStep";
    end
end