function [reward, ev, st]=getNextRew(rt,st)
reward=st.lookup(rt,st.sample(rt)+1);
ev=st.ev(rt);
st.sample(rt) = st.sample(rt) + 1;
