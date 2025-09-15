function figure_handle = visualize_vertical_forces(ttl, Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, Fv_INL, Fv_INR, t_FP, t_Mot)
    if (length(t_FP) == length(Fv_FP1)) && (length(t_Mot) == length(Fv_INL))
        figure_handle = figure;  
        p1 = plot(t_FP, Fv_FP1, 'red');
        t1 = "FP1";
        hold on
        p2 = plot(t_FP, Fv_FP2, 'green');
        t2 = "FP2";
        p3 = plot(t_FP, Fv_FP3, 'blue');
        t3 = "FP3";
        p4 = plot(t_FP, Fv_FP4, 'yellow');
        t4 = "FP4";
        p5 = plot(t_Mot, Fv_INL, 'magenta');
        t5 = "INL";
        p6 = plot(t_Mot, Fv_INR, 'cyan');
        t6 = "INR";
        legend([p1, p2, p3, p4, p5, p6], [t1, t2, t3, t4, t5, t6]);
        title(ttl)
        hold off
    else
        figure_handle = figure;  
        p1 = plot(Fv_FP1, 'red');
        t1 = "FP1";
        hold on
        p2 = plot(Fv_FP2, 'green');
        t2 = "FP2";
        p3 = plot(Fv_FP3, 'blue');
        t3 = "FP3";
        p4 = plot(Fv_FP4, 'yellow');
        t4 = "FP4";
        p5 = plot(Fv_INL, 'magenta');
        t5 = "INL";
        p6 = plot(Fv_INR, 'cyan');
        t6 = "INR";
        legend([p1, p2, p3, p4, p5, p6], [t1, t2, t3, t4, t5, t6]);
        title(ttl)
        hold off
    end
end