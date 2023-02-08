

template <>
inline int mpm::DamageMesh<1>::GetNodeRawIndex(Eigen::Matrix<int, 1, 1> index) {
    return index[1];
};
template <>
inline int mpm::DamageMesh<2>::GetNodeRawIndex(Eigen::Matrix<int, 2, 1> index) {
    return index[1] + (index[0] * mesh_size[1]);
};
template <>
inline int mpm::DamageMesh<3>::GetNodeRawIndex(Eigen::Matrix<int, 3, 1> index) {
    return index[2] + (index[1] * mesh_size[1]) + (index[0] * mesh_size[2] * mesh_size[1]);
};

template <unsigned Tdim>
template <typename Toper>
inline void mpm::DamageMesh<Tdim>::iterate_over_neighbours(Toper oper,mpm::ParticleBase<Tdim> & p,double distance){
    const int nodal_size = std::ceil(distance/resolution_);
    VectorDim index = PositionToIndex();
    if constexpr (Tdim == 1){
        for(int x = -nodal_size;x<=nodal_size;++x) {
            auto di = Eigen::Matrix<int,Tdim,1>() << x;
            auto & node = GetNode(index + di);
            node.iterate_over_particles(std::bind(oper,p,std::placeholders::_1));
        }
    }
    if constexpr (Tdim == 2){
        for(int x = -nodal_size;x<=nodal_size;++x) {
            for(int y = -nodal_size;y<=nodal_size;++y) {
                auto di = Eigen::Matrix<int,Tdim,1>() << x << y;
                auto & node = GetNode(index + di);
                node.iterate_over_particles(std::bind(oper,p,std::placeholders::_1));
            }
        }
    }
}
