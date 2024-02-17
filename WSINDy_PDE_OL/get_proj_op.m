function W = get_proj_op(W,UB)

    if ~isempty(UB)
  W = W.*(abs(W)<UB);
    end

end
