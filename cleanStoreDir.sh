rm prep/*
for i in {1..48}; do
	rm -r dir${i}/*
	rm store${i}/*
done
