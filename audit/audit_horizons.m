function report = audit_horizons(report)
N=numel(report.horizons);assert(N>=2,'v2:Horizon','Horizon audit requires >=2 horizons');
checks=[];
for j=1:N
    c=compare_horizon_results(report.baseline_cases{j},report.baseline_cases{end},[],[]);
    report.baseline_comparisons{j}=c;
    if j<N,checks(end+1)=c.screen_passed;end 
end
for k=1:numel(report.policy_cases)
    record=report.policy_cases{k};
    for j=1:N
        c=compare_horizon_results(record.results{j},record.results{end},report.baseline_cases{j},report.baseline_cases{end});
        record.comparisons{j}=c;
        if j<N,checks(end+1)=c.screen_passed;end 
    end
    report.policy_cases{k}=record;
end
report.horizon_screen_passed=all(checks);report.screen_failure_count=sum(~checks);
report.reference_horizon=report.horizons(end);report.infinite_horizon_accuracy_validated=false;

end
