function val = eval_counter(action, increment)
% EVAL_COUNTER
% Centralized evaluation counter for all limit-state calls

    persistent eval_count

    if isempty(eval_count)
        eval_count = 0;
    end

    switch action
        case 'reset'
            eval_count = 0;

        case 'add'
            eval_count = eval_count + increment;

        case 'get'
            % do nothing

        otherwise
            error('Unknown action');
    end

    val = eval_count;

end