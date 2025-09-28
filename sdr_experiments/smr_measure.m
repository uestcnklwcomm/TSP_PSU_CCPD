function SMR = smr_measure(x1,support1, x2, support2)

    
SMR = sum(abs(x1(support2)))/sum(abs(x1(support1))) + sum(abs(x2(support1)))/sum(abs(x2(support2)));
SMR = 0.5*SMR;

end