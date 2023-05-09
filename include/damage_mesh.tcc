

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
    const int nodal_size = std::ceil((2*distance)/resolution_);
    IndexDim index = PositionToIndex(p.coordinates());
    //In this case we are forced to do if constepxr or something much harder
    //We are not allowed to partially specialise template arguments, i.e. <Tdim=1,Toper=Any>
    //If constexpr still allows us to not check at runtime, but allows for this specialisation
    if constexpr (Tdim == 1){
        for(int x = -nodal_size;x<=nodal_size;++x) {
            IndexDim di = (IndexDim() << x).finished();
            auto & node = GetNode(index + di);
            //node.iterate_over_particles(std::bind(oper,&p,std::placeholders::_1));
        }
    }
    if constexpr (Tdim == 2){
        for(int x = -nodal_size;x<=nodal_size;++x) {
            for(int y = -nodal_size;y<=nodal_size;++y) {
                IndexDim di = (Eigen::Matrix<int, Tdim, 1>() << x, y).finished();
              if (InBounds(index + di)) {
                auto & node = GetNode(index + di);
                node.iterate_over_particles(std::bind(oper,&p,std::placeholders::_1));
              }
            }
        }
    }
    if constexpr (Tdim == 3){
        for(int x = -nodal_size;x<=nodal_size;++x) {
            for(int y = -nodal_size;y<=nodal_size;++y) {
            for (int z = -nodal_size; z <= nodal_size; ++z) {
              IndexDim di = (Eigen::Matrix<int, Tdim, 1>() << x, y, z).finished();
              if (InBounds(index + di)) {
                auto& node = GetNode(index + di);
                node.iterate_over_particles(
                    std::bind(oper, &p, std::placeholders::_1));
              }
            }
            }
        }
    }
}
