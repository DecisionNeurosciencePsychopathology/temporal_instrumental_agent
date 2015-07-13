function [reward, st]=getNextRew(rt,st)
reward=st.lookup(rt,st.sample(rt));
st.sample(rt) = st.sample(rt) + 1;