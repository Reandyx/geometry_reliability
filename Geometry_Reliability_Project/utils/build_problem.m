function problem = build_problem(a, b, model_type)

    switch model_type
        case 'global'
            problem = get_problem_synthetic(a, b);

        case 'local'
            problem = get_problem_local_curvature(a, b);

        otherwise
            error('Unknown model_type: %s', model_type);
    end

end