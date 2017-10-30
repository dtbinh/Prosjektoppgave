function checkNormalDistribution(figure_num, arr, arr_name)
    figure(figure_num);
    histogram(arr);
    hold on
    histogram(sqrt(var(arr))*(randn(1,length(arr))+mean(arr)));
    legend({arr_name, '$\mathcal{N}(\mu_{sig}, \sigma_{sig}^2)$'}, 'Interpreter', 'latex')
    hold off
    title('Histogram')
end