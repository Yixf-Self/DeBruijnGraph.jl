function de_bruijn_ize(str::ASCIIString, k::Integer)
    edges = Array{Tuple{ASCIIString,ASCIIString},1}()
    nodes = Set{ASCIIString}()
    for i in 1:(length(str)-k)
        push!(edges, (str[i:i+k-1],str[i+1:i+k]))
        push!(nodes, str[i:i+k-1], str[i+1:i+k])
    end
    return edges,nodes
end

type Node
    node_name::ASCIIString
    num_in::Int64
    num_out::Int64

    function Node(name::ASCIIString)
        new(name,0,0)
    end
end

isSemiBalanced(nd::Node) = abs(nd.num_in - nd.num_out) == 1
isBalanced(nd::Node) = nd.num_in == nd.num_out
name(nd::Node) = nd.node_name

pointer_from_node(nd::Node) = convert(Ptr{Node}, pointer_from_objref(nd))

function chop(str::ASCIIString, k::Int64)
    num_kmers = length(str)-k+1
    mers = Array{Tuple{ASCIIString,ASCIIString,ASCIIString},1}()
    for i in 1:num_kmers
        kmer = str[i:i+k-1]
        push!(mers,(kmer,kmer[1:end-1],kmer[2:end]))
    end
    mers
end

type DebruijnGraph
    G::Dict{Node,Array{Node,1}}
    nodes::Dict{ASCIIString,Node} # map k-1 mers to Node objects

    nsemi::Int64
    nbal::Int64
    nneither::Int64
    head::Ptr{Node}
    tail::Ptr{Node}
    
    function DebruijnGraph(strIter,k,circularize=false)
        G = Dict{Node,Array{Node,1}}()
        nodes = Dict{ASCIIString,Node}()
        for str in strIter
            if circularize
                str = string(str, str[1:k-1])
            end
            for (_,km1L,km1R) in chop(str,k)
                if !in(km1L, keys(nodes))
                    nodes[km1L] = Node(km1L)
                end
                if !in(km1R, keys(nodes))
                    nodes[km1R] = Node(km1R)
                end
                nodeL = nodes[km1L]
                nodeR = nodes[km1R]

                nodeL.num_in += 1
                nodeL.num_out += 1
                
                if !in(nodeL,keys(G))
                    G[nodeL] = Node[nodeR]
                else
                    push!(G[nodeL],nodeR)
                end
            end
        end
        nsemi,nbal,nneither = 0,0,0
        head,tail = Ptr{Node}(0),Ptr{Node}(0)
        
        for nd in values(nodes)
            if isBalanced(nd)
                nbal += 1
            elseif isSemiBalanced()
                if nd.num_in == nd.num_out + 1
                    tail = pointer_from_node(nd)
                else
                    head = pointer_from_node(nd)
                end
                nsemi += 1
            else
                nneither += 1
            end
        end
        
        new(G,nodes,nsemi,nbal,nneither,head,tail)        
    end

end

nnodes(g::DebruijnGraph) = length(g.nodes)
nedges(g::DebruijnGraph) = length(g.G)
hasEulerianWalk(g::DebruijnGraph) = g.nneither ==0 && g.nsemi == 2
hasEulerianCycle(g::DebruijnGraph) = g.nneither ==0 && g.nsemi == 0
isEulerian(g::DebruijnGraph) = hasEulerianWalk(g) || hasEulerianCycle(g)

function eulerianWalkOrCyle(g::DebruijnGraph)
    @assert isEulerian(g) == true
    
    G = g.G
    if hasEulerianWalk(g)
       if !in(g.head, G[unsafe_load(g.tail)])
           G[unsafe_load(g.tail)] = Node[unsafe_load(g.head)]
       else
           push!(G[unsafe_load(g.tail)],unsafe_load(g.head))
       end
    end

    tour = Node[]
    src = collect(keys(G))[1]

    function _visit(nd)
        while length(G[nd])>0
            dst = pop!(G[nd])
            _visit(dst)
        end
        push!(tour,nd)
    end

    _visit(src)
    tour =  reverse(tour)[1:end-1]

    if hasEulerianWalk(G)
        ind = findfirst(tour,unsafe_load(g.head))
        tour = vcat(tour[ind:end],tour[1:ind-1])
    end

    map(name,tour)
end
