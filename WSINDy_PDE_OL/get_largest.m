function y=get_largest(x,s)
    [~,I]=sort(abs(x),'descend');
    y=x*0;
    y(I(1:s))=x(I(1:s));
end