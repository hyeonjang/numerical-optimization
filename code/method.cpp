float Method::bisection(float start, float end)
{
    assert( function(start)*function(end)<0 );
    std::cout << "interval: [" << start << ", " << end << "]" << std::endl;
    std::cout << "function: [" << function(start) << ", " << function(end) << "]" << std::endl;

    auto midpoint = (start + end)/2.f;
    auto result = function(midpoint);

    std::cout << "midpoint: " << midpoint << std::endl;
    std::cout << "result: " << result << std::endl;

    if(result==0 || end-start<MIN)
        return midpoint;

    if(function(midpoint)*function(start)<0)
        midpoint = bisection(start, midpoint);
    else
        midpoint = bisection(midpoint, end);
    
    return midpoint;
}

